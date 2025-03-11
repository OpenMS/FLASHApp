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


![Configure Parameters](/static/Images/flashdeconv_upload.png)

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

 
![Configure Parameters](/static/Images/flashdeconv_configure.png)

---

## **3️⃣ Running the Workflow**

1. Click on the **"Run"** tab.
2. Set **log details** to `minimal` (or another level as needed).
3. Click **"Start Workflow"** to begin the deconvolution process.
4. Monitor the **log output** to track progress.

  
![Run Workflow](/static/Images/flashdeconv_run.png)

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


![Manual Upload](/static/Images/flashdeconv_manual_upload.png)

---

## **6️⃣ Using Example Data**

1. Click the **"Example Data"** tab.
2. Click **"Load Example Data"** to use the preloaded dataset.
3. The example data will appear in the uploaded experiments table.

 
![Example Data](/static/Images/flashdeconv_example_data.png)

---

## **7️⃣ Layout Manager**

The **Layout Manager** allows users to customize the experiment display settings.

1. Select the **number of experiments** to view at once.
2. Click **"Select..."** to choose components to add:
   - **MS1 raw heatmap**
   - **MS1 deconvolved heatmap**
   - **Scan table**
   - **Deconvolved spectrum**
   - **Annotated spectrum**
   - **Mass table**
   - **2D plot**
3. Click **Save** to apply changes.


![Layout Manager](/static/Images/flashdeconv_layout_manager.png)

---

## **8️⃣ Viewing Results in FLASHViewer**

1. Navigate to the **Viewer** tab.
2. Choose an experiment from the dropdown.
3. View the selected one from Layout manager: scan table, mass table, annotated spectrum, and deconvolved spectrum etc.

![FLASHViewer](/static/Images/flashdeconv_viewer.png)

---

## **9️⃣ Downloading Results**

1. Navigate to the **Download** tab.
2. Locate the experiment you want to download.
3. Click **"Prepare Download"** to generate the downloadable files.
4. To delete an experiment, click the **trash icon** next to the experiment name.

![Download Results](/static/Images/flashdeconv_download.png)

---



# **FLASHTnT Guide**

## **1️⃣ Uploading MS Data & Database**

1. Navigate to **FLASHTnT > Workflow** in the sidebar.
2. Click **"File Upload"** to upload your `.mzML` file.

![Download Results](/static/Images/flashTnT_upload.png)

3. Click the **"Database"** tab to upload the necessary **FASTA** database files.
4. Click **"Add Database"** to confirm the upload.


![Download Results](/static/Images/flashTnT_databaseupload.png)
---
## **2️⃣ Configuring Parameters**
1. Click the **"Configure"** tab.
2. Select your uploaded **mzML file**.
3. Choose the **FASTA database** file.
4. There are two sub-tabs for configuring parameters: **FLASHDeconv** and **FLASHTnT**.
5. Adjust FLASHTnT parameters:
Modify Ex Parameters (max_mod_mass, max_mod_count).Adjust **general settings** such as:
   - Threads
   - FDR settings
   - Ion types (`b`, `y`)
   - Modification mass limits
   - Tag lengths
   
5. Click **Save** to apply settings.


![Download Results](/static/Images/flashtnt_configure.png)
![Download Results](/static/Images/flashTnT_configure2.png)


---

## **3️⃣ Running the Workflow**

1. Click on the **"Run"** tab.
2. Click **"Start Workflow"** to begin.
3. Monitor the progress in the log output.

![Download Results](/static/Images/flashtnt_run.png)

---

## **4️⃣ Layout Manager**

1. Navigate to **Layout Manager**.
2. Select an experiment to view at once e.g, 1,2,3,4 or 5.
3. Now based on experiment number select the component to be added like **Protein table, Sequence view, Internal Fragment Map, Tag table, and Spectrum View**.

![Download Results](/static/Images/layoutmanager_tnt.png)

---

## **5️⃣ Viewer**

1. Choose the experiment.
2. View the selected one from Layout manager.

![Download Results](/static/Images/FlashTnT_Viewer.png)

---

## **6️⃣ Manual Result Upload & Example Data**

1. Click **"Manual Result Upload"** to upload manually processed data.
2. Click **"Example Data"** to load a sample dataset.

![Download Results](/static/Images/manual_result_upload.png)

---

## **7️⃣ Downloading Results**

1. Navigate to the **Download** tab.
2. Locate the experiment you want to download.
3. Click **"Prepare Download"** to generate the downloadable files.
4. To delete an experiment, click the **trash icon** next to the experiment name.
 
![Download Results](/static/Images/download_tnt.png)

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
  📄 [Read on arXiv](https://arxiv.org/abs/2411.13189)

---

🚀 **You're now ready to use FLASHDeconv & FLASHTnT!**







