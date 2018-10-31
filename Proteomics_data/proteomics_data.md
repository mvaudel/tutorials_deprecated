# Proteomics data

In this tutorial, navigate proteomics data at various stages of a quantitative proteomic analysis. Please note that there is a vast variety of acquisition modes resulting in differences in the data. Notably, one should distinguish (i) _targeted_ proteomics, where specific compounds are targeted by the mass spectrometer; and (ii) _untargeted_ acquisition, where the mass spectrometer screens in an unbiased fashion. Two acquisition modes can be distinguished in _untargeted_ proteomics: (i) _Data Dependent Acquisition_, DDA, where the mass spectrometer attempts at fragmenting one compound at a time by selecting them based on mass, and (ii) _Data Independent Acquisition_, DIA, where the mass spectrometer fragments all compounds within a given mass range. Here, we will focus on untargeted DDA proteomics data.


## Raw data

Raw proteomics data consist of multiple mass spectrometry runs, each corresponding to the acquisition of a sample. A sample can be spread across multiple runs, which is referred to as _sample fractionation_. The organization of samples and files depends on the experimental design, and it is important to take the experimental design into account when analysing data.


## Public proteomics data

In a global effort for scientific transparency, scientists worldwide share the proteomics data supporting their conclusion. The ProteomeXchange consortium is a worldwide coordinated effort for sharing proteomics data [(1)](#references).

![](http://www.proteomexchange.org/px_members.png)

:wrench: Navigate to [ProteomeXchange](http://www.proteomexchange.org), click on _Access Data_.

![](http://www.proteomexchange.org/access_data.png)

In this section, you can navigate the data provided by the community and look for data sets of interest.

:wrench: Navigate to [ProteomeXchange](http://www.proteomexchange.org), click on _Access Data_.


## Format and tools

Mass spectrometers produce raw data in various formats which are not necessarily open. These formats can be converted to mzML, the reference standard for mass spectrometry files [(2)](#references). The reference tool to process raw files is Proteowizzard [(3)](#references), it allows navigating and converting most mass spectrometry formats.

:wrench: Install [Proteowizzard](http://proteowizard.sourceforge.net/)
* Download 







## References

(1) [ProteomeXchange provides globally coordinated proteomics data submission and dissemination](https://www.ncbi.nlm.nih.gov/pubmed/24727771)

(2) [mzML-a community standard for mass spectrometry data](https://www.ncbi.nlm.nih.gov/pubmed/20716697)

(3) [A cross-platform toolkit for mass spectrometry and proteomics](https://www.ncbi.nlm.nih.gov/pubmed/23051804)

