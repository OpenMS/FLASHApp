# **Top-down analysis in three steps**

FLASHApp turns top-down MS data into proteoform masses, identifications and, if you like, a comparison of your samples. Every analysis follows the same three pages under **Top-down analysis** in the sidebar:

1. **Describe experiment**: your MS runs, a protein database to identify proteoforms, and what you want to learn, in the **Your analysis** box.
2. **Run analysis**: one click, then follow the progress.
3. **Viewer**: every run, and every comparison of samples, laid out for the analysis you chose.

The short videos below walk through each step with example data. No knowledge of the tools behind FLASHApp is needed to follow them.

@video guide-1-start.mp4 | Start an analysis: choose your runs and describe the experiment

---

## **Which analysis do I need?**

Three questions decide the analysis. You answer them in the **Your analysis** box on **Describe experiment**.

| Question | Answer | What FLASHApp does |
|---|---|---|
| Do you want to know **which proteins and proteoforms** are in the sample? | **Identify proteoforms** on | Searches every run against your protein database. Off: the masses are found, nothing is identified (intact-mass or MS1-only runs, or no database). |
| What kind of **sample** is it? (only with identification) | **Native or purified proteins** | Settings for a few proteins per sample: native MS, a purified protein, targeted runs. |
| | **Cell lysate or complex mixture** | Settings for many proteins separated by LC before MS/MS. |
| Do you want to **compare samples**? | **None** | Every run is looked at on its own. |
| | **Label-free: compare runs** | Each sample is its own run. Runs are linked by mass and retention time and compared between conditions. Needs at least two runs, ideally three or more per condition. |
| | **Isobaric labels (TMT, iTRAQ)** | The samples were labelled and pooled. The reporter ion of each label gives its sample's abundance. Several runs of one pool count as fractions. |

That gives six analyses, plus the sample type whenever proteoforms are identified:

| Analysis | Identify | Quantification |
|---|---|---|
| Deconvolution only | off | None |
| Deconvolution + label-free quantification | off | Label-free |
| Deconvolution + TMT/iTRAQ quantification | off | Isobaric labels |
| Deconvolution + identification (native or lysate) | on | None |
| Deconvolution + identification + label-free quantification | on | Label-free |
| Deconvolution + identification + TMT/iTRAQ quantification | on | Isobaric labels |

Under the box, a short summary (🧭) says in plain words what the next start will do. Read it before you go on.

**Not sure whether your runs are native or lysate?** Once you select the MS data, FLASHApp looks at the runs (how many MS1 spectra, how they were acquired) and says which settings fit. A button applies them.

---

## **1️⃣ Describe experiment**

1. Open **Top-down analysis › Describe experiment**.
2. In **Your analysis**, pick the **MS data**. **📂 Browse** beside it opens your files:
   - **Your files** lists them by folder. Tick the runs to analyse, or pick a folder to select its runs at once. Type a folder name in the *Folder* column to sort a run into it; FLASHApp suggests folders from the file names (for example `Control_R1`, `Control_R2`, `Treated_R1`), and a run's folder becomes its condition under **Samples and groups**.
   - **Data on your computer** takes new runs: drag them in or browse for them. `.mzML`, Thermo `.raw` and Bruker timsTOF `.d` are accepted; vendor files are converted for you.
   - **Cloud data** imports them from a link (an object store, OneDrive, Zenodo or Figshare).
   - **Import from PRIDE** reanalyses a public dataset: enter its accession (for example `PXD019247`) and press **Look up**. Its MS runs go to **MS data** and a FASTA it includes to **Protein database**; untick what you do not need, then **Import**. When a run is there both as mzML and as a RAW file, only the mzML is ticked.

   Without files of your own, a small example data set is used.
3. Switch **Identify proteoforms** on or off, pick the **Protein database** (a **FASTA** or **PEFF** file, added the same way with **📂 Browse**, whose **UniProt proteome** tab also finds the proteome of your organism by name or taxonomy ID and downloads it; FLASHApp adds the decoy sequences it needs) and the **Sample** type, and choose the **Quantification** (see the table above).
4. The numbered tabs below follow the steps of the analysis. You rarely need more than the first:
   - **1 · Samples and groups**: which sample belongs to which condition (see below).
   - **2 · Deconvolution**, **3 · Identification**, **4 · Fragment annotation**: the settings of each tool. The sample type already sets them; change them only if you know why. **Show advanced parameters** shows every setting.
   - **5 · Linking runs** (label-free only): how runs are matched by mass and retention time.

### Samples and groups

- **Label-free**: one row per run. Give each run its **condition** (for example *Control* and *Treated*) and, if they exist, its replicate and fraction. A sample sheet (design.tsv or SDRF) can be imported instead.
- **Isobaric labels**: choose the **Labels** (TMT 6-, 10-, 11-plex, TMTpro 16-, 18-plex, iTRAQ 4-, 8-plex). The table then lists one row per label channel: give each channel its **sample** name and its **condition**. Every channel needs a condition to be compared; a channel you leave empty is shown but not tested.

---

## **2️⃣ Run analysis**

1. Open **Top-down analysis › Run analysis**. The box at the top repeats the analysis you chose.
2. Press **Start Workflow**. The log shows each run as it is processed. You can close the page and come back with the same workspace.
3. When it finishes, open the **Viewer**.

---

## **3️⃣ Viewer**

The Viewer shows one run at a time (more with **Edit layout**). Pick it in **choose experiment**. Everything on the page is linked: select a spectrum, a protein or a proteoform in one table and every spectrum, table and plot follows. **💡 Add tool results** shows the output of a FLASHDeconv or FLASHTnT run made outside FLASHApp.

@video guide-2-viewer.mp4 | Explore a run in the Viewer and arrange your own layout

Until you save your own layout, each run opens in the layout made for its analysis:

| Run | What you see first |
|---|---|
| Deconvolved only | The masses over retention time (heatmap), the spectra with their masses, the raw and deconvolved spectrum and a 3D signal-to-noise view. |
| Identified, native or purified proteins | The identified proteins, the sequence with its fragments and modifications, the 3D structure, and the sequence tags behind each identification. |
| Identified, lysate | Proteins, their proteoforms and the spectra matching each (PrSMs) side by side, the proteoform map over retention time and mass, the sequence view, fragment errors and quality control of the identifications. |
| Labelled with TMT or iTRAQ | The layout of its kind as above, with the **reporter ions** of every spectrum and of the selected spectrum, right below the table that selects the spectrum. |

**Edit layout** opens the Layout Manager: drag components into the grid, choose how many runs to show at once, and **Save**. **Reset** goes back to the defaults above. A saved layout is used for every run of the workspace.

### Components

Search the list in the Layout Manager by name, or expand its groups:

- **Heatmaps**: *MS1 raw heatmap* (m/z over retention time) and *MS1 deconvolved heatmap* (mass over retention time).
- **Tables**: *Scan table* (every spectrum), *Mass table* (the masses of the selected spectrum), *Protein table* (identified proteins, with a **Characterization** badge: *Characterized* when both ends are known, every modification is placed on one residue and at least 25% of the bonds are covered by fragments; *Ambiguous* below 10%; *Partial* in between), *Proteoforms*, *PrSMs*, *Tag table*.
- **Spectra & plots**: *Raw spectrum*, *Deconvolved spectrum*, *Isotope envelope fit* (observed against expected isotope pattern: a shifted top points to a wrong monoisotopic mass), *Precursor XIC by charge*, *3D S/N plot*, *Run QC*.
- **Sequence**: *Sequence view* (fragments and modifications on the sequence, internal fragments as a purple strip), *Structure view (3D)* (the AlphaFold model, or a structure you upload, with the cleavage sites), *Internal fragment map*, *Fragment mass errors*.
- **Quantification**: *Reporter ions per spectrum (TMT, iTRAQ)* and *Reporter ions of the selected spectrum*.

The **Sequence view** of a run that was not identified shows a sequence you set in the view itself with **Change sequence** (the DNA icon above the sequence); it is saved for the workspace.

---

## **5️⃣ Comparing samples**

With a quantification chosen, the run ends with a comparison of the samples. It appears in **choose experiment** as **Comparison: …** (the number of runs, or the labels and samples, and the date). Selecting it shows a summary, the comparison table, a volcano plot, the abundance profile of the selected entry, every sample's value with its detection status, a PCA of the samples, a heatmap and quality control.

**Edit comparison** changes the conditions, the level compared (features, proteoforms or proteins; the latter two need identification) and the statistics (fold-change and adjusted p-value cutoffs, filtering, imputation, normalization, test). The comparison is recomputed in seconds; nothing is run again.

@video guide-3-label-free.mp4 | Compare conditions without labels: one run per sample

@video guide-4-tmt.mp4 | Compare TMT-labelled samples: reporter ions and the comparison

**How the numbers come about**

- **Label-free**: the abundance of a mass is the area of its signal in each run. Runs are matched by mass and retention time; where a mass was not detected, its signal there is read (match between runs). With identification, each linked mass gets the proteoform the run identified.
- **Isobaric labels**: the abundance in a sample is the reporter ion intensity of its channel in the fragmentation spectra, summed over the spectra of the same mass (and over the fractions). With identification, spectra are summed per proteoform and protein (at 1% false discovery rate).

---

## **6️⃣ Download results**

**Download** lists every run. **Prepare Download** packs its result files (deconvolved spectra, identifications, quantification tables); the trash icon deletes a run from the workspace.

---

## **Words used in FLASHApp**

- **Deconvolution**: turning the many charge states and isotope peaks of a protein in a spectrum into one neutral mass.
- **Proteoform**: one exact molecular form of a protein, with its modifications and cleavages.
- **PrSM**: a proteoform-spectrum match, one spectrum explained by one proteoform.
- **q-value / FDR**: the expected share of wrong identifications. 1% is the usual cutoff.
- **Sequence tag**: a short run of amino acids read directly from the fragment masses, used to find the protein in the database.
- **Reporter ion**: the small ion each isobaric label releases on fragmentation; its intensity is the sample's abundance.
- **Feature**: a mass followed over retention time (label-free), or the spectra of one mass (labelled).
