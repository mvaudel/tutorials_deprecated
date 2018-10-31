# Proteomics data

In this tutorial, navigate proteomics data at various stages of a quantitative proteomic analysis. Please note that there is a vast variety of acquisition modes resulting in differences in the data. Notably, one should distinguish (i) _targeted_ proteomics, where specific compounds are targeted by the mass spectrometer; and (ii) _untargeted_ acquisition, where the mass spectrometer screens in an unbiased fashion. Two acquisition modes can be distinguished in _untargeted_ proteomics: (i) _Data Dependent Acquisition_, DDA, where the mass spectrometer attempts at fragmenting one compound at a time by selecting them based on mass, and (ii) _Data Independent Acquisition_, DIA, where the mass spectrometer fragments all compounds within a given mass range. Here, we will focus on untargeted DDA proteomics data.


## Raw data

Raw proteomics data consist of multiple mass spectrometry runs, each corresponding to the acquisition of a sample. A sample can be spread across multiple runs, which is referred to as _sample fractionation_. The organization of samples and files depends on the experimental design, and it is important to take the experimental design into account when analysing data.


## Public proteomics data

In a global effort for scientific transparency, scientists worldwide share the proteomics data supporting their conclusion. The ProteomeXchange consortium is a worldwide coordinated effort for sharing proteomics data [(1)](#references).

![ProteomeXchange overview](http://www.proteomexchange.org/px_members.png "ProteomeXchange overview")

:pencil2: Navigate to [ProteomeXchange](http://www.proteomexchange.org), click on `Access Data`.

In this section, you can navigate the data provided by the community and look for data sets of interest.

:pencil2: Open the data set `PXD001819`.

You now see information on the data set, the publication as well as metadata on the project.

:pencil2: Click `PRIDE project URI`.

You are now redirected to PRIDE [(1)](#references), the repository where the data were deposited. You can access more information as well as the original files.

:pencil2: Click `Download Project Files`.

You now see the files as they were provided by the submitter.

:pencil2: Scroll down to `RAW Files` and download `UPS1_12500amol_R1.raw`.


## Format and tools

Mass spectrometers produce raw data in various formats which are not necessarily open. These formats can be converted to mzML, the reference standard for mass spectrometry files [(3)](#references). The reference tool to process raw files is Proteowizzard [(4)](#references), it allows navigating and converting most mass spectrometry formats. Note that the libraries needed to open the files are often provided for Windows only, Mac and Linux users will encounter difficulties working with the raw data.

:pencil2: Install [Proteowizzard](http://proteowizard.sourceforge.net/)


## Navigating raw data

Navigating the raw data will allow you to find specific spectra, and look at the spectra acquired by the instrument. It is the first step of quality control (QC) and can help spotting many mistakes.

:pencil2: Go to the installation folder and open `SeeMS`.

:pencil2: Open the file `UPS1_12500amol_R1.raw` downloaded in the previous section.

You should see the following picture.

![SeeMS Overview](Proteomics_data/images/seeMS_1.jpg?raw=true "SeeMS Overview").








## References

(1) [ProteomeXchange provides globally coordinated proteomics data submission and dissemination](https://www.ncbi.nlm.nih.gov/pubmed/24727771)

(2) [PRIDE: the proteomics identifications database](https://www.ncbi.nlm.nih.gov/pubmed/16041671)

(3) [mzML-a community standard for mass spectrometry data](https://www.ncbi.nlm.nih.gov/pubmed/20716697)

(4) [A cross-platform toolkit for mass spectrometry and proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23051804)

