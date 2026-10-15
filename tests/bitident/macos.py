"""Bit-identity checks for macOS builds (CI job "bit-identity (macos ...)").

The macOS counterpart of tests/bitident/windows.py and of the Linux haswell job.  Run at the
repository root after the workflow has extracted the pre-constexpr header (BASELINE) to old/.
For each build (compiler and flags for this machine's architecture) it compiles
tests/bitident/dump.cpp at -O2 against the current header and against old/, and checks:

  * old vs new, raw bits: must be identical in every build;
  * non-contracting builds, canonical hashes vs the glibc goldens in golden.inc, with
    "<type>@macos-<arch>" entries taking precedence where present (abs/normalize call the
    C library's hypot) and long double compared as double where it is binary64 (arm64);
  * every build vs the first (plain) build of the same compiler, listing which operations
    change; builds with a "same" entry must be identical to that build.

Compilers: Apple Clang (xcrun clang++, the SDK's libc++) and, with --llvm PREFIX, Homebrew
LLVM clang with its own libc++.  Exit status 1 on any failure; Markdown summary to --report.
"""
import argparse, os, platform, re, subprocess, sys

COMMON = ["-std=c++23", "-O2"]

def builds(arch, apple, llvm):
    # (name, compiler, extra flags, library code contracted?, must equal this build)
    b = []
    # the toolchain's clang++ called directly needs the SDK (xcrun or CMake would add it)
    sdk = subprocess.run(["xcrun", "--show-sdk-path"], capture_output=True, text=True).stdout.strip()
    comps = [("apple-clang", apple, ["-isysroot", sdk] if sdk else [])]
    if llvm:
        comps.append(("llvm-clang", llvm[0], llvm[1]))
    for cname, cc, ld in comps:
        if arch == "arm64":
            # AArch64: fma is baseline and Clang contracts within an expression by default
            b += [(cname, cc, ld, True, None),
                  (cname + " -ffp-contract=off", cc, ld + ["-ffp-contract=off"], False, None)]
        else:
            b += [(cname, cc, ld, False, None),
                  (cname + " -march=haswell", cc, ld + ["-march=haswell"], True, None),
                  (cname + " -march=haswell -ffp-contract=off", cc, ld + ["-march=haswell", "-ffp-contract=off"], False, cname)]
    return b

def load_golden(path):
    g = {}
    for m in re.finditer(r'\{"([^"]+)", "([^"]+)", 0x([0-9a-f]+)ULL\}', open(path).read()):
        g[(m.group(1), m.group(2))] = m.group(3)
    return g

def parse(text):
    d = {}
    for line in text.splitlines():
        t, op, h = line.split()
        d[(t, op)] = h
    return d

def compile_run(cc, flags, inc, out, mode):
    cmd = [cc] + COMMON + ["-I" + inc, "tests/bitident/dump.cpp", "-o", out] + flags
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return None, (r.stdout + r.stderr)[-3000:]
    r = subprocess.run([out] + ([mode] if mode else []), capture_output=True, text=True)
    if r.returncode != 0:
        return None, "run failed: %d %s" % (r.returncode, r.stderr[-2000:])
    return parse(r.stdout), ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apple", default="clang++")
    ap.add_argument("--llvm", default="", help="Homebrew LLVM prefix")
    ap.add_argument("--report", default="bitident-macos.md")
    a = ap.parse_args()
    arch = "arm64" if platform.machine() in ("arm64", "aarch64") else "x86_64"
    llvm = None
    if a.llvm and os.path.exists(a.llvm + "/bin/clang++"):
        p = a.llvm
        llvm = (p + "/bin/clang++", ["-L" + p + "/lib/c++", "-L" + p + "/lib/unwind", "-lunwind", "-Wl,-rpath," + p + "/lib/c++"])
    golden = load_golden("tests/bitident/golden.inc")
    os.makedirs("bi", exist_ok=True)
    fails, lines = [], ["## Bit identity (macos, %s)" % arch, ""]
    if not llvm:
        lines.append("Homebrew LLVM not found (%s): llvm-clang builds skipped.\n" % (a.llvm or "-"))
    canon = {}
    lines += ["| Build | hashes | old vs new (raw) | vs Linux goldens (canonical) | vs plain build of the same compiler |", "|---|---|---|---|---|"]
    for i, (name, cc, flags, contracts, same) in enumerate(builds(arch, a.apple, llvm)):
        tag = "bi/b%d" % i
        new_raw, err = compile_run(cc, flags, "include/quaternion", tag + "n", "")
        if new_raw is None:
            fails.append("%s: build/run failed" % name); lines.append("| `%s` | build failed | | | |" % name)
            print("::group::%s error\n%s\n::endgroup::" % (name, err)); continue
        new_c, _ = compile_run(cc, flags, "include/quaternion", tag + "c", "canonical")
        canon[name] = new_c
        old_raw, err = compile_run(cc, flags, "old/quaternion", tag + "o", "")
        if old_raw is None:
            ovn = "old header does not build"
            print("::group::%s old header error\n%s\n::endgroup::" % (name, err))
            fails.append("%s: old header does not build" % name)
        else:
            diff = sorted(k for k in new_raw if old_raw.get(k) != new_raw[k])
            ovn = "identical" if not diff else "DIFFER: " + ", ".join("%s %s" % k for k in diff)
            if diff: fails.append("%s: old vs new differ: %s" % (name, diff))
        if contracts:
            vg = "skipped (contracting build)"
        else:
            bad = []
            ld64 = arch == "arm64"  # long double is binary64 on Apple arm64
            for (t, op), h in new_c.items():
                tt = "double" if (t == "long_double" and ld64) else t
                cand = [tt + "@macos-" + arch, tt]
                want = next((golden[(c, op)] for c in cand if (c, op) in golden), None)
                if want != h: bad.append("%s %s %s (golden %s)" % (t, op, h, want))
            vg = "all %d match" % len(new_c) if not bad else "%d DIFFER: %s" % (len(bad), "; ".join(bad))
            if bad: fails.append("%s: golden mismatch: %s" % (name, bad))
        plain = same or name.split(" ")[0]
        if plain in canon and plain != name:
            d = sorted("%s %s" % k for k in new_c if canon[plain].get(k) != new_c[k])
            vp = ("identical" if not d else "%d differ: %s" % (len(d), ", ".join(d))) + " (vs `%s`)" % plain
            if same and d: fails.append("%s: expected identical to %s: %s" % (name, plain, d))
        else:
            vp = "(plain build)"
        lines.append("| `%s` | %d | %s | %s | %s |" % (name, len(new_raw), ovn, vg, vp))
    lines += ["", "Failures: %d" % len(fails)] + ["* " + f for f in fails]
    open(a.report, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(main())
