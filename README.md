# mmSYGNAL-risk-prediction-models
ML risk prediction models generated from mmSYGNAL

# **mmSYGNAL risk prediction**  
  
Tutorial for applying mmSYGNAL risk models to a new patient diagnosed with 
multiple myeloma. Risk models were generated separately based on the genetic
abnormality subtypes patients in our training cohort exhibited. The agnostic 
model was generated with all patients in training cohort regardless of subtype
and is appropriate for patients that do not exhibit a relevant genetic subtype
or when subtype information is not obtainable.  
  
### **mmSYGNAL risk models**    
  
See analysis/risk_model_tutorial.nb.html for code.  
   
The following risk models are available.    
  
no subtype    
   
1. agnostic  

cytogentic subtypes  
  
2. amp(1q)  
3. del(13)  
4. del(1p) 
5. t(4;14)
  
RNAseq subtypes  
  
5. FGFR3  
  
### **Steps required for pipeline**    
  
1. Generate gene expression data for patient and annotate with ensembl gene ids. 
RNAseq data is recommended but mmSYGNAL also functions with microarray data.
  
2. Apply mmSYGNAL model to gene expression data to generate program activity. 
See 'Wall et al, Genetic program activity delineates risk, relapse, and therapy responsiveness in multiple myeloma,
Precision Oncology, 2021' for instructions on how to apply mmSYGNAL.
  
3. Format patient's program activity into a data.frame where the row is the
patient's program activity column titles are the program label (0 to 140).  

4. Identify which genetic abnormality risk subtypes the patient exhibits and
apply the appropriate risk models. See code below for examples. 
  
5. If desired choose only the highest grade models (A > B > C) if a patient 
exhbits multiple subtypes. If a patient shows multiple subtypes of the same
grade then take the mean of the risk scores.  
  
6. The patients risk classification is generated based on their predicted 
risk score as shown below:  
  
+ low risk: < 0.5  
  
+ high risk: >= 0.5 and < 0.6  
  
+ extreme risk: >= 0.6  
  
