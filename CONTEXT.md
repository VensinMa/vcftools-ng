# vcftools-ng compatibility engine

This context defines the project language used to preserve VCFtools 0.1.17
behavior while sharing variant decoding across parallel analyses.

## Language

**Compatibility oracle**:
The original VCFtools 0.1.17 executable whose complete output bytes define
correct behavior.
_Avoid_: Reference implementation, expected approximation

**Real subset**:
The fixed 2,300,000-record fixture containing 100,000 records from each of 23
chromosomes and all 412 samples.
_Avoid_: Toy fixture, full dataset

**Golden**:
A complete output artifact produced by the compatibility oracle from the
BGZF VCF encoding of the real subset; exact BCF container tests use the BCF
encoding when the legacy VCF-to-BCF path cannot represent the real header.
_Avoid_: Snapshot sample, expected rows

**Accepted variant**:
A variant that has passed every site, genotype, sample, and ordered thinning
decision for the active invocation.
_Avoid_: Decoded record, candidate site

**Analysis payload**:
Batch-owned, column-oriented GT/DP and contribution data decoded once for use
by every requested analysis.
_Avoid_: Per-output cache, SiteResult vectors

**Ordered commit**:
The single input-order transition that applies thinning, marks accepted
variants, and feeds deterministic writers and reducers.
_Avoid_: Worker merge, unordered reduction

**Input backend**:
The container-specific mechanism that turns plain VCF, indexed BGZF/BCF, or
unindexed BGZF into the same ordered shard stream.
_Avoid_: File extension branch, decompression thread count

**Ordered shard**:
A uniquely owned, stable-ordinal span of complete records that may be read and
decoded out of order but is exposed to exact consumers in original order.
_Avoid_: Arbitrary chunk, unordered region

**Record-start ownership**:
The rule assigning an indexed record to the one half-open coordinate window
containing its start position, even if index overlap semantics return it from
adjacent windows.
_Avoid_: Deduplicate after output, assume index queries never overlap

**Deterministic analysis lane**:
A reducer that permanently owns an output axis and observes accepted variants
in original site order so floating-point accumulation matches the oracle.
_Avoid_: Arbitrary worker partial, nondeterministic reduction

**Exact compatibility**:
Complete-file byte equality with a golden for both 8-thread and 16-thread
executions.
_Avoid_: Numerically close, semantically equivalent
