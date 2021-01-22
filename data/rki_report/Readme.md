# RKI Report-file Elements

[toc]

## Source -Please check first

* https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/DESH.html
* https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Cryptshare-Anleitung.pdf?__blob=publicationFile [v1.1 (2021-01-21)]

## Data Overview
**The report must contain following seven elements in this order:**

1. IMS_ID
2. SENDING_LAB
3. DATE_DRAW
4. SEQ_TYPE
5. SEQ_REASON
6. SAMPLE_TYPE
7. OWN_FASTA_ID

In the following a short description of each element is given.


### 1. IMS_ID

Sequencing-based pseudonym as unique identifier for the aggregation in DEMIS ("Deutschen Elektronischen Melde- und Informationssystem für den Infektionsschutz"). 

**Format: IMS-12345-CVDP-00001**
* IMS: permanent prefix 
* 12345: 5-digit identifier of the sequencing laboratory ("Untersuchungslabor"), analog to the already existing DEMIS-system (DEMIS-10001 to currently DEMIS-10563). The list is managed by the DEMIS-Geschäftsstelle. If you are not registered already pls reach out to demis@rki.de
* CVDP: 4-digit DEMIS-abbreviation, which is directly dedicated to the "Meldetatbestand". Later additional pathogen-abbreviations, aside from SARS-CoV-2, applied.
* 00001: Ongoing Number, which in phase 0 is autonmously continued by the laboratory (later in phase 1 the system will automatically gernerate it).


### 2. SENDING_LAB

12345: 5-digit identifier of the sending laboratory, analog to the already existing DEMIS-system.

**ATTENTION:**
This only applys for laboratories, which don´t sequence on their own, but instead sent their samples to other laboratories for sequencing. In case that the sending lab is also the sequencing lab the digit from the IMS_ID and the SENDING_LAB-id can be identical.


### 3. DATE_DRAW

Date of the sample isolation n ISO8601 (YYYYMMDD)


### 4. SEQ_TYPE

Used sequencing-platform. "OXFORD_NANOPORE" is provided automatically as entry.


### 5. SEQ_REASON

Cause for the sequencing. Choose one entry from the following list:

|Entry|Description|
|-|-|
|X|Unknown to the sequencing laboratory|
|N|No (e. g. random selection of a PCR-positive sample for sequencing)|
|Y|Yes, but the kind of mutation or variante is unknown (to the sequencing laboratory) 
|A|Yes, it exists evidence for the mutation/variante from previous diagnostic [spezifying in textfield after entry-letter]|

**Note for "A":**
* Textfield, max. length 64 signs
* Mutation to specify in "[ ]", in case of multiple mutations divided by "/"
* Example entry: A[B.1.1.7/B.1.351]


### 6. SAMPLE_TYPE

Type of sample. Choose one entry from the following list:

|Entry|Description|
|-|-|
|s001|Upper respiratory swab sample (specimen)|
|s002|Nasopharyngeal swab (specimen)|
|s003|Swab from nasal sinus (specimen)|
|s004|Anterior nares swab (specimen)|
|s005|Oropharyngeal aspirate (specimen)|
|s006|Nasopharyngeal aspirate (specimen)|
|s007|Lower respiratory sample (specimen)|
|s008|Bronchoalveolar lavage fluid sample (specimen)|
|s009|Sputum specimen (specimen)|
|s010|Specimen from trachea obtained by aspiration (specimen)|
|s011|Pleural fluid specimen (specimen)|
|s012|Specimen from lung obtained by biopsy (specimen)
|s013|Blood specimen (specimen)|
|s014|Plasma specimen or serum specimen or whole blood specimen (specimen)|
|s015|Whole blood sample (specimen)|
|s016|Stool specimen (specimen)|
|s017|Urine specimen (specimen)|
|s018|Lower respiratory fluid sample (specimen)|
|s019|Nasopharyngeal washings (specimen)|
|s020|Plasma specimen (specimen)|
|s021|Saliva specimen (specimen)|
|s022|Serum specimen (specimen)|
|s023|Specimen unsatisfactory for evaluation (finding)|
|s024|Swab of internal nose (specimen)|
|s025|Throat swab (specimen)|
|X|Unknown (to the sequencing laboratory)

Value Set is geared to SNOMED CT and SNOMED CT COVID-19 Related Content(https://simplifier.net/covid-19labormeldung/materialsarscov2)

### 7. OWN_FASTA_ID

Laboratory-internal identifier, which enables the distinct assignment of the FASTA-file to the sequence (given in the FASTA-header). Autoprovided by PoreCov-workflow, according to your sample-names.

**NOTE:**
Used in phase 0 for the assignment of the metadata from the .csv-file to the sequence-data in the FASTA-file.