"""Bit-identity checks for Windows builds (CI job "bit-identity (windows ...)").

Run from a Visual Studio developer shell at the repository root.  For each build in BUILDS
(compiler and flags for this --arch) it compiles tests/bitident/dump.cpp twice, against the
current header and against the pre-constexpr header at BASELINE (extracted to old/ by the
workflow), and runs both.  It then checks:

  * old vs new, raw bits: must be identical in every build (as in the Linux haswell job);
  * non-contracting builds, canonical hashes vs the glibc goldens in golden.inc, with the
    Universal CRT's own entries for abs/normalize ("double@ucrt"; "double@ucrt-arm64" on
    ARM64 if present) and long double compared as double;
  * every build vs the first (plain) build of the same compiler, listing which operations
    change; builds with a "same" entry must be identical to that build.

MSVC with /fp:contract: since the headers compile their own code with fp_contract(off) there
(see quaternion.hpp), the library is not contracted and the new header must give the plain
build's bits and the goldens.  The pre-constexpr header has no such pragma, so old vs new
is expected to differ in that build (only the transcendental-free operations that the
compiler contracts); that difference is listed, not failed.

Exit status 1 on any failure.  A Markdown summary goes to --report.
"""
import argparse, os, re, shutil, subprocess, sys

COMMON = ["/nologo", "/std:c++latest", "/EHsc", "/utf-8", "/O2"]

def builds(arch, clangcl):
    # (name, compiler, extra flags, library code contracted?, must equal this build)
    b = []
    if arch == "x64":
        b += [("msvc", "cl", [], False, None),
              ("msvc /arch:AVX2", "cl", ["/arch:AVX2"], False, "msvc"),
              ("msvc /arch:AVX2 /fp:contract", "cl", ["/arch:AVX2", "/fp:contract"], False, "msvc")]
        if clangcl:
            b += [("clang-cl", clangcl, [], False, None),
                  ("clang-cl /arch:AVX2", clangcl, ["/arch:AVX2"], True, None),
                  ("clang-cl /arch:AVX2 -ffp-contract=off", clangcl, ["/arch:AVX2", "/clang:-ffp-contract=off"], False, "clang-cl")]
    else:
        b += [("msvc", "cl", [], False, None),
              ("msvc /fp:contract", "cl", ["/fp:contract"], False, "msvc")]
        if clangcl:
            b += [("clang-cl", clangcl, [], True, None),
                  ("clang-cl -ffp-contract=off", clangcl, ["/clang:-ffp-contract=off"], False, None)]
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
    exe = out + ".exe"
    cmd = [cc] + COMMON + flags + ["/I" + inc, "tests/bitident/dump.cpp", "/Fe:" + exe, "/Fo:" + out + ".obj"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return None, (r.stdout + r.stderr)[-3000:]
    r = subprocess.run([exe] + ([mode] if mode else []), capture_output=True, text=True)
    if r.returncode != 0:
        return None, "run failed: %d %s" % (r.returncode, r.stderr[-2000:])
    return parse(r.stdout), ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", choices=["x64", "arm64"], required=True)
    ap.add_argument("--clang-cl", default="")
    ap.add_argument("--report", default="bitident-windows.md")
    a = ap.parse_args()
    clangcl = a.clang_cl if a.clang_cl and os.path.exists(a.clang_cl) else ""
    golden = load_golden("tests/bitident/golden.inc")
    os.makedirs("bi", exist_ok=True)
    fails, lines = [], ["## Bit identity (windows, %s)" % a.arch, ""]
    if not clangcl:
        lines.append("clang-cl not found (%s): clang-cl builds skipped.\n" % (a.clang_cl or "-"))
    canon = {}
    lines += ["| Build | hashes | old vs new (raw) | vs Linux goldens (canonical) | vs plain build of the same compiler |", "|---|---|---|---|---|"]
    for i, (name, cc, flags, contracts, same) in enumerate(builds(a.arch, clangcl)):
        tag = "bi/b%d" % i
        new_raw, err = compile_run(cc, flags, "include/quaternion", tag + "n", "")
        if new_raw is None:
            fails.append("%s: build/run failed" % name); lines.append("| %s | build failed | | | |" % name)
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
            if diff and "/fp:contract" in name:
                ovn = "differ as expected (old header contracted; new: fp_contract off): " + ", ".join("%s %s" % k for k in diff)
            elif diff:
                fails.append("%s: old vs new differ: %s" % (name, diff))
        if contracts:
            vg = "skipped (contracting build)"
        else:
            bad = []
            for (t, op), h in new_c.items():
                tt = "double" if t == "long_double" else t  # long double is binary64 here
                cand = [tt]
                if tt == "double" and op in ("abs", "normalize"):
                    cand = ["double@ucrt-arm64", "double@ucrt"] if a.arch == "arm64" else ["double@ucrt"]
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
