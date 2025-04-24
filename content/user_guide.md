# **FLASHDeconv & FLASHTnT User Guide**

Welcome to the **FLASHDeconv & FLASHTnT User Guide**. This guide provides a step-by-step walkthrough on using **FLASHDeconv** and **FLASHTnT**, from uploading data to processing and viewing results.

---

## **1️⃣ Uploading MS Data**

1. Navigate to **FLASHDeconv > Workflow** in the sidebar.
2. Click **"File Upload"** to upload your `.mzML` file.
3. **Two options to add files:**
   - **Drag and Drop** your file into the upload box.
   - **Browse Files** to manually select the `.mzML` file.
4. Click **"Add MS Data"** to confirm the upload.

![Upload](/static/images/flashdeconv_upload.png)


---

## **2️⃣ Configuring Parameters**

1. Click the **"Configure"** tab.
2. Select your uploaded **mzML file** (should appear in the file list).
3. Adjust **parameters**:
   - Set **threads** (recommended: 8).
   - Choose **General Settings** (e.g., `keep_empty_out`, `write_detail`).
   - Configure **FD settings** like `report_FDR`, `merging_method`, `min_mass`, `max_mass`, `min_charge`.
   - Adjust **SD settings** for deconvolution accuracy.
4. If you want to know more about each parameter go to this link: https://openms.de/FLASHDeconv

 
![Configure Parameters](/static/images/flashdeconv_configure.png)

---

## **3️⃣ Running the Workflow**

1. Click on the **"Run"** tab.
2. Set **log details** to `minimal` (or another level as needed).
3. Click **"Start Workflow"** to begin the deconvolution process.
4. Monitor the **log output** to track progress.

  
![Run Workflow](/static/images/flashdeconv_run.png)

---

## **4️⃣ Viewing Results**

1. Once the workflow is finished, check the **log messages**.
2. Navigate to the **"Viewer"** tab in the sidebar to analyze the deconvoluted data.
3. If needed, **download results** by clicking **"Download Files"**.
   (Will be explained later in step 7 and 8 in this guide)

---

## **5️⃣ Manual Result Upload**

1. Click on the **"Manual Result Upload"** tab.
2. Upload FLASHDeconv output files (`*_annotated.mzML` & `*_deconv.mzML`) or TSV files for ECDF Plot.
3. Browse files or **drag and drop** them into the upload section.
4. Click **"Add files to workspace"** to finalize.


![Manual Upload](/static/images/flashdeconv_manual_upload.png)

---

## **6️⃣ Using Example Data**

1. Click the **"Example Data"** tab.
2. Click **"Load Example Data"** to use the preloaded dataset.
3. The example data will appear in the uploaded experiments table.

 
![Example Data](/static/images/flashdeconv_example_data.png)

---

## **7️⃣ Layout Manager**

The **Layout Manager** allows users to customize the experiment display settings.

1. Select the **number of experiments** to view at once.
2. Click **"Select..."** to choose components to add:
   - **MS1 raw heatmap** - 2D heatmap of raw MS signals with m/z (y-axis), retention time (x-axis), and intensity as color gradient.
   - **MS1 deconvolved heatmap** - Displays raw MS signals as a 2D heatmap with monoisotopic mass (y-axis), retention time (x-axis), and intensity as a color gradient.
   - **Scan table** - Lists scan details (e.g., number, retention time, precursor mass). 
   - **Deconvolved spectrum** - Plots deconvolved spectrum for a scan (summed intensity vs. monoisotopic mass).
   - **Raw spectrum** - Plots raw spectrum for a scan(intensity vs m/z).
   - **Mass table** - Display deconvolved masses for a selected scan with properties.
   - **3D S/N plot** - Visualizes S/N ratio of deconvolved masses in 3D.
3. Click **Save** to apply changes.


![Layout Manager](/static/images/flashdeconv_layout_manager.png)

---

## **8️⃣ Viewing Results in FLASHViewer**

1. Navigate to the **Viewer** tab.
2. Choose an experiment from the dropdown.
3. View the selected one from Layout manager: scan table, mass table, annotated spectrum, and deconvolved spectrum etc.

![FLASHViewer](/static/images/flashdeconv_viewer.png)

---

## **9️⃣ Downloading Results**

1. Navigate to the **Download** tab.
2. Locate the experiment you want to download.
3. Click **"Prepare Download"** to generate the downloadable files.
4. To delete an experiment, click the **trash icon** next to the experiment name.

![Download Results](/static/images/flashdeconv_download.png)

---



# **FLASHTnT Guide**

## **1️⃣ Uploading MS Data & Database**

1. Navigate to **FLASHTnT > Workflow** in the sidebar.
2. Click **"File Upload"** to upload your `.mzML` file.

![Download Results](/static/images/flashTnT_upload.png)

3. Click the **"Database"** tab to upload the necessary **FASTA** database files.
4. Click **"Add Database"** to confirm the upload.


![Download Results](/static/images/flashTnT_databaseupload.png)
---
## **2️⃣ Configuring Parameters**
1. Click the **"Configure"** tab.
2. Select your uploaded **mzML file**.
3. Choose the **FASTA database** file.
4. There are two sub-tabs for configuring parameters: **FLASHDeconv** and **FLASHTnT**.
5. Adjust FLASHTnT parameters:
Adjust **general settings** such as:
   - **FDR settings** (prsm_fdr, pro_fdr) – Set thresholds for precursor and proteoform-level FDR (e.g., 1.00%).
   - **Ion types** (ion_type) – Ion series to consider when generating tags (e.g., b, y).
   - **Max modification count** (max_mod_count) – Maximum number of allowed modifications per protein.

   
5. Click **Save** to apply settings.


![Download Results](/static/images/flashtnt_configure.png)
![Download Results](/static/images/flashTnT_configure2.png)


---

## **3️⃣ Running the Workflow**

1. Click on the **"Run"** tab.
2. Click **"Start Workflow"** to begin.
3. Monitor the progress in the log output.

![Download Results](/static/images/flashtnt_run.png)

---

## **4️⃣ Layout Manager**

1. Navigate to Layout Manager.
2. Select how many experiments to display simultaneously (e.g., 1–5).
3. Choose **components** to add for each experiment:
- **Protein Table** – Lists proteins identified by FLASHTnT, including accession, modifications, and score.
- **Sequence View** – Annotates tags, PTMs, and fragments. Visualizes FLASHTnT results; for FLASHDeconv data must be user-supplied.
- **Internal Fragment Map** – Shows internal fragment ions from the selected scan
- **Tag Table** – Lists sequence tags with corresponding information.
- **Spectrum View** – Shows the annotated spectrum with matched peaks.

![Download Results](/static/images/layoutmanager_tnt.png)

---

## **5️⃣ Viewer**

1. Choose the experiment.
2. View the selected one from Layout manager.

![Download Results](/static/images/FlashTnT_Viewer.png)

---

## **6️⃣ Manual Result Upload & Example Data**

1. Click **"Manual Result Upload"** to upload manually processed data.
2. Click **"Example Data"** to load a sample dataset.

![Download Results](/static/images/manual_result_upload.png)

---

## **7️⃣ Downloading Results**

1. Navigate to the **Download** tab.
2. Locate the experiment you want to download.
3. Click **"Prepare Download"** to generate the downloadable files.
4. To delete an experiment, click the **trash icon** next to the experiment name.
 
![Download Results](/static/images/download_tnt.png)

---


## **📖 Need Help?**

If you have any questions or need assistance, feel free to contact our support team.

### **FLASHApp Support Contacts:**
- **Tom Müller**: [tom.mueller@uni-tuebingen.de](mailto:tom.mueller@uni-tuebingen.de)
- **Ayesha Feroz**: [ayesha.feroz@uni-tuebingen.de](mailto:ayesha.feroz@uni-tuebingen.de)

---

## **📚 Relevant Publications**
For more information about the research behind FLASHDeconv & FLASHTnT, refer to the following publications:

- **Jeong, K., Kim, J., Gaikwad, M., Hidayah, S. N., Heikaus, L., Schlüter, H., & Kohlbacher, O.** (2020).  
  *FLASHDeconv: Ultrafast, high-quality feature deconvolution for top-down proteomics.*  
  **Cell Systems, 10(2), 213-218.e6**  
  📄 [Read the paper](https://doi.org/10.1016/j.cels.2020.01.003)

- **Müller, T. D., Siraj, A., Walter, A., Kim, J., Wein, S., von Kleist, J., Feroz, A., Pilz, M., Jeong, K., Sing, J. C., Charkow, J., Röst, H. L., & Sachsenberg, T.** (2024, November 20).  
  *OpenMS webapps: Building user-friendly solutions for MS analysis.*  
  📄 [Read the paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00872)

- **Müller, T., Kim, J., Almaguer, A., et al.** (2025, April 17).
 *FLASHApp: Interactive Data Analysis and Visualization for Top-Down Proteomics.*
 📄 [Read the paper](https://www.authorea.com/users/914240/articles/1287405-flashapp-interactive-data-analysis-and-visualization-for-top-down-proteomics?commit=1447908dbdd26d9a2312890c9c400d96f2b171f7)
---

🚀 **You're now ready to use FLASHDeconv & FLASHTnT!**







