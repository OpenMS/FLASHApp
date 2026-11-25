# Processing report: flashtnt

|  |  |
|---|---|
| App | FLASHApp 2.0.0-alpha.11 |
| Workflow | flashtnt |
| Started | 2026-11-14 20:06 |
| pyOpenMS | 3.6.0 |
| Python | 3.11.18 |
| Report written | 2026-11-14 20:07 |

## Software

| Tool | Version |
|---|---|
| FLASHDeconv | unknown |
| DecoyDatabase | unknown |
| FLASHTnT | unknown |
| FLASHAnnotator | unknown |
| FLASHQuantConsensus | 3.6.0-pre-FDmerge-2026-11-09 |

## Settings

| Setting | Value |
|---|---|
| identify | true |
| quantification | label-free |
| few_proteins | false |
| fasta-file | example_database.fasta |

## Input files

Every tool that worked on each input file, in order. A tool call that took several files at once is listed under each of them. The parameters of each call are listed under **Tool calls**.

### Label-free_SpikeIn-1x_R1.mzML

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 1 | FLASHDeconv | unknown | done | 2 s | out.tsv, spec1.tsv, spec2.tsv, spec3.tsv, spec4.tsv, out_deconv.mzML, quant.tsv, anno_annotated.mzML, toppic_ms1.msalign, toppic_ms2.msalign, toppic_ms1.feature, toppic_ms2.feature |
| 3 | FLASHTnT (with 1 other file) | unknown | done | 6 s | tags.tsv, protein.tsv, prsms.tsv |
| 4 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |
| 21 | FLASHQuantConsensus (with 5 other files) | 3.6.0-pre-FDmerge-2026-11-09 | done | 0 s | consensus.tsv, identifications.tsv |

### example_database.fasta

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 2 | DecoyDatabase | unknown | done | 0 s | database.fasta |
| 3 | FLASHTnT (with 1 other file) | unknown | done | 6 s | tags.tsv, protein.tsv, prsms.tsv |
| 4 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |
| 7 | DecoyDatabase | unknown | done | 0 s | database.fasta |
| 8 | FLASHTnT (with 1 other file) | unknown | done | 7 s | tags.tsv, protein.tsv, prsms.tsv |
| 9 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |
| 12 | DecoyDatabase | unknown | done | 0 s | database.fasta |
| 13 | FLASHTnT (with 1 other file) | unknown | done | 6 s | tags.tsv, protein.tsv, prsms.tsv |
| 14 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |
| 17 | DecoyDatabase | unknown | done | 0 s | database.fasta |
| 18 | FLASHTnT (with 1 other file) | unknown | done | 8 s | tags.tsv, protein.tsv, prsms.tsv |
| 19 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |

### Label-free_SpikeIn-1x_R2.mzML

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 6 | FLASHDeconv | unknown | done | 3 s | out.tsv, spec1.tsv, spec2.tsv, spec3.tsv, spec4.tsv, out_deconv.mzML, quant.tsv, anno_annotated.mzML, toppic_ms1.msalign, toppic_ms2.msalign, toppic_ms1.feature, toppic_ms2.feature |
| 8 | FLASHTnT (with 1 other file) | unknown | done | 7 s | tags.tsv, protein.tsv, prsms.tsv |
| 9 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |
| 21 | FLASHQuantConsensus (with 5 other files) | 3.6.0-pre-FDmerge-2026-11-09 | done | 0 s | consensus.tsv, identifications.tsv |

### Label-free_SpikeIn-0.1x_R1.mzML

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 11 | FLASHDeconv | unknown | done | 2 s | out.tsv, spec1.tsv, spec2.tsv, spec3.tsv, spec4.tsv, out_deconv.mzML, quant.tsv, anno_annotated.mzML, toppic_ms1.msalign, toppic_ms2.msalign, toppic_ms1.feature, toppic_ms2.feature |
| 13 | FLASHTnT (with 1 other file) | unknown | done | 6 s | tags.tsv, protein.tsv, prsms.tsv |
| 14 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |
| 21 | FLASHQuantConsensus (with 5 other files) | 3.6.0-pre-FDmerge-2026-11-09 | done | 0 s | consensus.tsv, identifications.tsv |

### Label-free_SpikeIn-0.1x_R2.mzML

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 16 | FLASHDeconv | unknown | done | 2 s | out.tsv, spec1.tsv, spec2.tsv, spec3.tsv, spec4.tsv, out_deconv.mzML, quant.tsv, anno_annotated.mzML, toppic_ms1.msalign, toppic_ms2.msalign, toppic_ms1.feature, toppic_ms2.feature |
| 18 | FLASHTnT (with 1 other file) | unknown | done | 8 s | tags.tsv, protein.tsv, prsms.tsv |
| 19 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |
| 21 | FLASHQuantConsensus (with 5 other files) | 3.6.0-pre-FDmerge-2026-11-09 | done | 0 s | consensus.tsv, identifications.tsv |

### out.tsv

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 21 | FLASHQuantConsensus (with 5 other files) | 3.6.0-pre-FDmerge-2026-11-09 | done | 0 s | consensus.tsv, identifications.tsv |

### prsms.tsv

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 21 | FLASHQuantConsensus (with 5 other files) | 3.6.0-pre-FDmerge-2026-11-09 | done | 0 s | consensus.tsv, identifications.tsv |

## Tool calls

Parameters that differ from the tool's defaults.

### Step 1: FLASHDeconv

Version: unknown · Status: done · Time: 2 s

| Parameter | Value |
|---|---|
| iq:type | none |
| quant | true |

### Step 2: DecoyDatabase

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| method | shuffle |
| shuffle_decoy_ratio | 1 |
| enzyme | no cleavage |

### Step 3: FLASHTnT

Version: unknown · Status: done · Time: 6 s

| Parameter | Value |
|---|---|
| native | false |

### Step 4: FLASHAnnotator

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| algorithm:internal_ions | false |
| ion_type | b, y |

### Step 5: Processing results

Version: unknown · Status: done · Time: 1 s

Parameters were not recorded for this call.

### Step 6: FLASHDeconv

Version: unknown · Status: done · Time: 3 s

| Parameter | Value |
|---|---|
| iq:type | none |
| quant | true |

### Step 7: DecoyDatabase

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| method | shuffle |
| shuffle_decoy_ratio | 1 |
| enzyme | no cleavage |

### Step 8: FLASHTnT

Version: unknown · Status: done · Time: 7 s

| Parameter | Value |
|---|---|
| native | false |

### Step 9: FLASHAnnotator

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| algorithm:internal_ions | false |
| ion_type | b, y |

### Step 10: Processing results

Version: unknown · Status: done · Time: 1 s

Parameters were not recorded for this call.

### Step 11: FLASHDeconv

Version: unknown · Status: done · Time: 2 s

| Parameter | Value |
|---|---|
| iq:type | none |
| quant | true |

### Step 12: DecoyDatabase

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| method | shuffle |
| shuffle_decoy_ratio | 1 |
| enzyme | no cleavage |

### Step 13: FLASHTnT

Version: unknown · Status: done · Time: 6 s

| Parameter | Value |
|---|---|
| native | false |

### Step 14: FLASHAnnotator

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| algorithm:internal_ions | false |
| ion_type | b, y |

### Step 15: Processing results

Version: unknown · Status: done · Time: 0 s

Parameters were not recorded for this call.

### Step 16: FLASHDeconv

Version: unknown · Status: done · Time: 2 s

| Parameter | Value |
|---|---|
| iq:type | none |
| quant | true |

### Step 17: DecoyDatabase

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| method | shuffle |
| shuffle_decoy_ratio | 1 |
| enzyme | no cleavage |

### Step 18: FLASHTnT

Version: unknown · Status: done · Time: 8 s

| Parameter | Value |
|---|---|
| native | false |

### Step 19: FLASHAnnotator

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| algorithm:internal_ions | false |
| ion_type | b, y |

### Step 20: Processing results

Version: unknown · Status: done · Time: 0 s

Parameters were not recorded for this call.

### Step 21: FLASHQuantConsensus

Version: 3.6.0-pre-FDmerge-2026-11-09 · Status: done · Time: 0 s

All parameters at their defaults.
