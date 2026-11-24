# Processing report: flashtnt

|  |  |
|---|---|
| App | FLASHApp 2.0.0-alpha.11 |
| Workflow | flashtnt |
| Started | 2026-11-14 20:06 |
| pyOpenMS | 3.6.0 |
| Python | 3.11.18 |
| Report written | 2026-11-14 20:06 |

## Software

| Tool | Version |
|---|---|
| FLASHDeconv | unknown |
| DecoyDatabase | unknown |
| FLASHTnT | unknown |
| FLASHAnnotator | unknown |

## Settings

| Setting | Value |
|---|---|
| identify | true |
| quantification | none |
| few_proteins | true |
| fasta-file | example_database.fasta |

## Input files

Every tool that worked on each input file, in order. A tool call that took several files at once is listed under each of them. The parameters of each call are listed under **Tool calls**.

### Native-identified_GapA_HCD.mzML

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 1 | FLASHDeconv | unknown | done | 1 s | out.tsv, spec1.tsv, spec2.tsv, spec3.tsv, spec4.tsv, out_deconv.mzML, quant.tsv, anno_annotated.mzML, toppic_ms1.msalign, toppic_ms2.msalign, toppic_ms1.feature, toppic_ms2.feature |
| 3 | FLASHTnT (with 1 other file) | unknown | done | 20 s | tags.tsv, protein.tsv, prsms.tsv |
| 4 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |

### example_database.fasta

| Step | Tool | Version | Status | Time | Wrote |
|---|---|---|---|---|---|
| 2 | DecoyDatabase | unknown | done | 0 s | database.fasta |
| 3 | FLASHTnT (with 1 other file) | unknown | done | 20 s | tags.tsv, protein.tsv, prsms.tsv |
| 4 | FLASHAnnotator (with 1 other file) | unknown | done | 0 s | annotated.tsv, fragments.tsv |

## Tool calls

Parameters that differ from the tool's defaults.

### Step 1: FLASHDeconv

Version: unknown · Status: done · Time: 1 s

| Parameter | Value |
|---|---|
| iq:type | none |

### Step 2: DecoyDatabase

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| method | shuffle |
| shuffle_decoy_ratio | 10 |
| enzyme | no cleavage |

### Step 3: FLASHTnT

Version: unknown · Status: done · Time: 20 s

| Parameter | Value |
|---|---|
| native | true |

### Step 4: FLASHAnnotator

Version: unknown · Status: done · Time: 0 s

| Parameter | Value |
|---|---|
| algorithm:internal_ions | true |
| ion_type | b, y |

### Step 5: Processing results

Version: unknown · Status: done · Time: 0 s

Parameters were not recorded for this call.
