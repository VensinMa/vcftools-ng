# Release workflow

This workflow keeps vcftools-ng releases reproducible without forcing every
release to have identical prose or benchmarks. Expensive release qualification
is separate from the permanent three-scenario development gate.

## 1. Freeze a candidate

- Choose the release version and a short highlight subtitle.
- Finish the three-scenario 1/2/4/8/16/32-thread development gate.
- Update the binary/CMake version and build the exact candidate that will be
  benchmarked and packaged.
- Record the source commit or candidate tree identity. Do not change runtime
  code after the final matrix without rerunning affected gates.

## 2. Correctness and safety

- Run the full CTest suite.
- Run the full ASan/UBSan suite.
- Run any new real-data golden gates introduced by the release.
- Preserve new Original goldens locally and commit their SHA-256 identities,
  commands, timings, and compact evidence.
- Document intentionally inherited Original defects and combinations rejected
  because Original output is corrupt or undefined.

## 3. Full seven-scenario qualification

Run only after explicit release authorization:

1. BGZF VCF + valid TBI/CSI;
2. BGZF VCF + automatic CSI;
3. BGZF VCF + `--no-auto-index`;
4. Plain VCF;
5. BCF + valid CSI;
6. BCF + automatic CSI;
7. BCF + `--no-auto-index`.

The standard host matrix uses 1/2/4/8/16/32 threads, strict serial execution,
and five repeats. Original is measured per distinct input format; scenarios
that differ only in vcftools-ng index policy share that format's Original
baseline. Automatic-index rows start from an independent path without a
sidecar on every run, so their time includes CSI construction.

The first candidate repeat for every scenario/thread must pass byte comparison
before repeats 2–5 begin. The release driver may additionally compare every
repeat. Inputs, original indexes, actual golden outputs, and output hashes must
remain available after the run. Large artifacts stay local; compact TSV,
manifest, environment, and reproduction scripts are committed.

When the release owner explicitly chooses staged qualification, run the driver
with `GATE_ONLY=1`. Publication may proceed only after all 42 first-repeat
candidate gates pass. Documentation and Release notes must then label every
value as single-run, state that repeats 2–5 are pending, and avoid five-run
mean or monotonic-scaling claims. Resume later with `GATE_ONLY=0`; completed
Original and first-repeat records are hash-validated and skipped. Commit the
five-run summary to master when it completes.

## 4. User-facing records

- Update `README.md` and `README.zh-CN.md` together.
- Keep visible language-switch links at the top of both files.
- Update the parameter compatibility matrix, version history, technical
  version record, release notes, benchmark summary, and known limitations.
- Include commands, test workload, Original policy, CPU/RSS, speedups,
  exactness, host environment, and links to machine-readable results.

## 5. Portable package

- Build the Linux x86_64 archive on the CentOS 7/manylinux2014 baseline.
- Verify its SHA-256.
- Test extraction, `--version`, bundled bcftools, runtime dependencies, a
  real-data compatibility smoke test, and automatic CSI in clean CentOS 7 and
  Ubuntu 20.04 containers.
- Use `packaging/linux-x86_64/test-portable.sh` so both clean-container checks
  stay identical across releases.
- Keep the archive layout (`bin`, `lib`, `libexec`) intact.

## 6. Publish

1. Commit all code, bilingual documentation, scripts, and compact evidence to
   direct `master`.
2. Push and verify `origin/master`.
3. Create the annotated release tag.
4. Create a GitHub Release named
   `vX.Y.Z — Highlight subtitle`.
5. Upload the portable archive and checksum.
6. Verify the public release body, assets, download links, source commit,
   README language links, and displayed version.

For a staged qualification release, the tag records the fully gated runtime
candidate while the later repeat-only evidence is a documentation update on
master. Runtime code must not change during the continuation; otherwise the
affected release matrix must restart from a newly frozen candidate.

If any validation or public verification fails, fix it and repeat the affected
stage. A release is complete only when source, documentation, benchmark
evidence, package, tag, and GitHub Release all describe the same version.
