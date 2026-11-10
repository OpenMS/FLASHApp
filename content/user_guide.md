# **Top-down analysis in three steps**

FLASHApp turns top-down MS data into proteoform masses, identifications and, if you like, a comparison of your samples. Every analysis follows the same three pages under **Top-down analysis** in the sidebar:

1. **Describe experiment**: your MS runs, a protein database to identify proteoforms, and what you want to learn, in the **Your analysis** box.
2. **Run analysis**: press **Run analysis** at the end of *Describe experiment*; the run starts right away and you follow its progress.
3. **Explore results**: every run, and every comparison of samples, laid out for the analysis you chose.

**Export results** (**Download results** in the web app) takes the result files with you, and **Match fragments** under **Other tools** matches fragment masses against a sequence without any run.

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
2. In **Your analysis**, pick the **MS data**. **Select a folder** above it selects every run of a folder at once. **📂 Browse** beside it opens your files:
   - **Your files** lists them by folder. Tick the runs to analyse, or pick a folder to select its runs at once. Type a folder name in the *Folder* column to sort a run into it; FLASHApp suggests folders from the file names (for example `Control_R1`, `Control_R2`, `Treated_R1`), and a run's folder becomes its condition under **Samples and groups**.
   - **Data on your computer** takes new runs: drag them in or browse for them. `.mzML`, Thermo `.raw` and Bruker timsTOF `.d` are accepted; vendor files are converted for you.
   - **Cloud data** imports them from a link (an object store, OneDrive, Zenodo or Figshare).
   - **Import from PRIDE** reanalyses a public dataset: enter its accession (for example `PXD055510`, intact TMT-labelled proteins) and press **Look up**. Its MS runs go to **MS data** and a FASTA it includes to **Protein database**; untick what you do not need, then **Import**. When a run is there both as mzML and as a RAW file, only the mzML is ticked.

   Without files of your own, a small example data set is used. The **Explore results** already holds the results of every example, each named after the analysis it shows: *Example: deconvolution only* (an LC-MS/MS run of a protein standard mix), *Example: native, identified* (a native top-down spectrum of E. coli GapA), *Example: lysate, identified* and *Example: label-free comparison (spike-in 1x vs 0.1x)* (four short runs of a six-protein mix spiked into an E. coli lysate at two amounts, from the FLASHQuant SpikeIn dataset, MassIVE MSV000091923), and *Example: TMT 6-plex, deconvolution only* and *Example: TMT comparison (8x vs 4x vs 1x)* (five minutes of one run of an E. coli lysate labelled at its cysteines with iodoTMT sixplex and mixed 8:4:1:1:4:8, PRIDE PXD055510, Matzanke, Kaulich et al., J Proteome Res 2025).
3. Switch **Identify proteoforms** on or off, pick the **Protein database** (a **FASTA** or **PEFF** file, added the same way with **📂 Browse**, whose **UniProt proteome** tab also finds the proteome of your organism by name or taxonomy ID and downloads it; FLASHApp adds the decoy sequences it needs) and the **Sample** type, and choose the **Quantification** (see the table above). Identification is optional: LC-MS/MS (DDA) runs can be deconvolved only, or compared label-free, without a database. The **Sample** type is always offered; it sets how proteoforms are identified, so without identification it changes nothing until **Identify proteoforms** is switched on.
4. The numbered tabs below follow the steps of the analysis you chose; a step the analysis does not run has no tab:
   - **Samples and groups** (with quantification only): which sample belongs to which condition (see below).
   - **Deconvolution**, and with identification **Identification** and **Fragment annotation**: the settings of each tool. The sample type already sets them; change them only if you know why. **All N parameters** beside each tool shows every setting.
   - **Linking runs** (label-free only): how runs are matched by mass and retention time.
5. Name the start in the field beside **Run analysis** (it is named after its files when left empty). The name is what **Explore results** and **Export results** list it as.

### Samples and groups

- **Label-free**: one row per run. Give each run its **condition** (for example *Control* and *Treated*) and, if they exist, its replicate and fraction. A sample sheet (design.tsv or SDRF) can be imported instead.
- **Isobaric labels**: choose the **Labels** (TMT 6-, 10-, 11-plex, TMTpro 16-, 18-plex, iTRAQ 4-, 8-plex). The table then lists one row per label channel: give each channel its **sample** name and its **condition**. Every channel needs a condition to be compared; a channel you leave empty is shown but not tested.

---

## **2️⃣ Run analysis**

1. Press **Run analysis** at the end of **Describe experiment**. The **Run analysis** page opens and the analysis starts at once; the box at the top repeats the analysis you chose. (**Start Workflow** on that page starts it again later.)
2. The log shows each run as it is processed. You can close the page and come back with the same workspace.
3. When it finishes, press **Explore results**.

---

## **3️⃣ Explore results**

The row at the top of **Explore results** says what is shown:

- **choose experiment** lists every start: the runs, and a comparison as **Comparison: …**. Each example is listed under its own name.
- **Single run | Comparison** (for a comparison only) shows one of its runs, or the comparison of all of them (see *Comparing samples* below). For a run that belongs to no comparison, **Comparison** offers to process more files to compare it with.
- **History** shows how the result came about: its settings and its processing report.
- **Add analysis** processes more files with the settings shown, or with **Add tool results** adds the output of a FLASHDeconv or FLASHTnT run made outside FLASHApp.

Everything on the page is linked: select a spectrum, a protein or a proteoform in one table and every spectrum, table and plot follows.

@video guide-2-viewer.mp4 | Explore a run and arrange your own view

### Views

The bar below this row switches between views. **Deconvolution**, **Native** and **Lysate** are ready-made; until you pick one, each run opens in the view made for its analysis:

| Run | What you see first |
|---|---|
| Deconvolved only (**Deconvolution**) | The masses over retention time (heatmap), the spectra with their masses, the raw and deconvolved spectrum and a 3D signal-to-noise view. |
| Identified, native or purified proteins (**Native**) | The identified proteins, the sequence with its fragments and modifications, the 3D structure, and the sequence tags behind each identification. |
| Identified, lysate (**Lysate**) | Proteins, their proteoforms and the spectra matching each (PrSMs) side by side, the proteoform map over retention time and mass, the sequence view, fragment errors, quality control of the identifications and the GO enrichment of the proteins you select. |
| Labelled with TMT or iTRAQ | The view of its kind as above, with the **reporter ions** of every spectrum and of the selected spectrum, right below the table that selects the spectrum. |

A view that needs results a run does not have (identification, for example) is greyed out; clicking it offers to run the missing analysis. A selection stays when you switch views.

**+** at the end of the bar makes a view of your own, starting as a copy of the current one: drag components into the grid, choose how many runs to show side by side, and save. **✏️** orders the views, keeps some out of the bar (they move to **More**) and edits or deletes your own. The ready-made views are never changed: editing one makes a copy.

### Components

Search the list in the view editor (**+** or **✏️**) by name, or expand its groups:

- **Heatmaps**: *MS1 raw heatmap* (m/z over retention time) and *MS1 deconvolved heatmap* (mass over retention time).
- **Tables**: *Scan table* (every spectrum), *Mass table* (the masses of the selected spectrum), *Protein table* (identified proteins, with a **Characterization** badge: *Characterized* when both ends are known, every modification is placed on one residue and at least 25% of the bonds are covered by fragments; *Ambiguous* below 10%; *Partial* in between), *Proteoforms*, *PrSMs*, *Tag table*.
- **Spectra & plots**: *Raw spectrum*, *Deconvolved spectrum*, *Isotope envelope fit* (observed against expected isotope pattern: a shifted top points to a wrong monoisotopic mass), *Precursor XIC by charge*, *3D S/N plot*, *Run QC*.
- **Sequence**: *Sequence view* (fragments and modifications on the sequence, internal fragments as a purple strip), *Structure view (3D)* (the AlphaFold model, or a structure you upload, with the cleavage sites), *Internal fragment map*, *Fragment mass errors*.
- **Quantification**: *Reporter ions per spectrum (TMT, iTRAQ)* and *Reporter ions of the selected spectrum*.
- **Functional analysis**: *GO enrichment of the selected proteins* (see below).

### GO enrichment of selected proteins

The GO enrichment component tests which Gene Ontology terms (biological process, cellular component, molecular function) are over-represented in the proteins you select, against all proteins identified in the run. Select several at once:

- in the *Proteins*, *Proteoforms*, *PrSMs* or *Protein table*, turn on the check boxes with the check mark button of the table's toolbar and tick rows (the box in the header ticks every row shown; Shift-click ticks a range),
- or pick box or lasso select in the mode bar of the *Proteoform heatmap* and drag over it. Shift adds to the selection.

Every component marks the selected proteins, and **Clear selection** (in the GO component, or the eraser in a plot's mode bar) starts over. The bar chart shows the most significant terms; how many is set next to the category. The GO terms of each protein are looked up once at MyGene.info by its UniProt accession (this needs a network connection) and kept in the workspace. The test is a one-sided Fisher's exact test, adjusted for the number of terms (Benjamini-Hochberg).

The **Sequence view** of a run that was not identified shows a sequence you set in the view itself with **Change sequence** (the DNA icon above the sequence); it is saved for the workspace.

---

## **Comparing samples**

With a quantification chosen, the run ends with a comparison of the samples. It appears in **choose experiment** as **Comparison: …** (the number of runs, or the labels and samples, and the date) and opens on **Comparison**; **Single run** looks at one of its runs instead. The comparison shows a summary, the comparison table, a volcano plot, the abundance profile of the selected entry, every sample's value with its detection status, a PCA of the samples, a heatmap, quality control and the GO enrichment of the entries you select.

The bar above the comparison switches between four views (and any of your own, made with **+** or ✏️):

- **Identification**: how many proteoforms (features, for runs not identified) each run detected itself and how many were matched from the other runs, in how many runs each was found, and what every two runs share. Below: every proteoform with its status run by run. Click one to list its PrSMs in every run, and click a PrSM to see its sequence view and annotated spectrum (from the run's own results; **Spectrum: yes** marks the PrSMs whose run is in this workspace). It needs no conditions.
- **Quantification** (the default): what changed between the conditions, as described above.
- **QC: label-free**: the consensus features: where they lie by retention time and mass (coloured by how many runs detected them; click one), the intensity of the values each run detected and of those matched into it, every feature with its charges and how many runs detected, matched or miss it, and the picked feature run by run. Then what each run detected and matched, how well the matches between runs scored, the total intensity of each run, the intensity distributions before and after normalization, how the runs correlate, the variation between replicates (CV), missing values, imputation and p-values.
- **QC: TMT, iTRAQ**: the reporter ions: their intensity per label (with the share of spectra that carry it), their ratio to the first label in the same spectrum, how many labels each spectrum carries, and every labelled spectrum with its reporter ions (pick the run when the set was split into several). Then the total intensity of each channel (unequal loading shows here), the distributions before and after normalization, how the channels correlate and cluster (PCA), the variation between replicates, missing values, imputation and p-values.

To test GO terms, select entries in the volcano plot (press **All changed**, **Higher in …** or **Lower in …**, or use box or lasso select in its mode bar), tick rows in the comparison table after turning its check boxes on, or select rows on the clustered heatmap with box or lasso select. The background is every protein of the entries tested.

**Edit comparison** changes the conditions, which runs are compared (the **Compared** column), the level compared (features, proteoforms or proteins; the latter two need identification) and the statistics (fold-change and adjusted p-value cutoffs, filtering, imputation, normalization, test). The comparison is recomputed in seconds; nothing is run again.

@video guide-3-label-free.mp4 | Compare conditions without labels: one run per sample

@video guide-4-tmt.mp4 | Compare TMT-labelled samples: reporter ions and the comparison

**How the numbers come about**

- **Label-free**: the abundance of a mass is the area of its signal in each run. Runs are matched by mass and retention time; where a mass was not detected, its signal there is read (match between runs). With identification, each linked mass gets the proteoform the run identified.
- **Isobaric labels**: the abundance in a sample is the reporter ion intensity of its channel in the fragmentation spectra, summed over the spectra of the same mass (and over the fractions). With identification, spectra are summed per proteoform and protein (at 1% false discovery rate).

---

## **4️⃣ Export results**

**Export results** lists every start. A run holds the deconvolved masses and spectra (TSV, mzML, msalign and feature files), with **Identify proteoforms** on also the proteins, PrSMs, tags and annotated fragments, and the parameters it ran with.

- In the desktop app, **Export** copies a run's files to the **Export folder** (one folder per run), and **Export all N runs** copies every run.
- In the web app, the page is called **Download results**: **Download** gives one ZIP per run, or every run at once.

✏️ renames a run, **History** shows how it came about and **Files** lists its files; the trash icon deletes it from the workspace. The run shown on **Explore results** can also be exported (or downloaded) there.

---

## **🧩 Match fragments**

**Other tools › Match fragments** matches a list of fragment masses against a sequence, without a deconvolution run:

1. Under **Observed masses**, say whether the masses are monoisotopic or average and neutral (M) or protonated (MH+), then **Paste** them (mass and intensity per line), load them **From a file** (text, tsv, csv or .msalign) or take them **From a run** of the workspace. The precursor mass is optional. **Load an example** fills in a scan of the lysate example.
2. Press **Match these masses** (**Match the masses of this file**, **Match this scan**).
3. Set the sequence in the sequence view with **Change sequence** (plain or ProForma), and add modifications by right-clicking a residue. The matched fragments are shown on the sequence.

---

## **Words used in FLASHApp**

- **Deconvolution**: turning the many charge states and isotope peaks of a protein in a spectrum into one neutral mass.
- **Proteoform**: one exact molecular form of a protein, with its modifications and cleavages.
- **PrSM**: a proteoform-spectrum match, one spectrum explained by one proteoform.
- **q-value / FDR**: the expected share of wrong identifications. 1% is the usual cutoff.
- **Sequence tag**: a short run of amino acids read directly from the fragment masses, used to find the protein in the database.
- **Reporter ion**: the small ion each isobaric label releases on fragmentation; its intensity is the sample's abundance.
- **Feature**: a mass followed over retention time (label-free), or the spectra of one mass (labelled).
