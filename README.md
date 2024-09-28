# PrecICE

## A unified workflow for data-driven precision cell fate engineering via highly multiplexed gene control


### Installation

**Using pip:**

Ensure you have Python and pip installed on your system. Then run the following command:

```pip install -r model/requirements.txt```

**Using conda:**

Ensure you have Anaconda or Miniconda installed on your system. Create a new conda environment (optional but recommended):

```
conda create --name cell_reprogram_env
conda activate cell_reprogram_env
```

Install the dependencies:
```
conda install --file model/requirements.txt
```





### Sample workflow


Worfklow for processing new dataset and running the model

#### Step 1. Initialize PRESCIENT data processing workflow and preprocess dataset
```
workflow = prescient_data(adata=adata,
                          dir='./workflow_dir', 
                          path='./hesc_dataset.h5ad')
```

#### Step 2. Computing differential expression
```
## Add cell type transition information
workflow.set_transition(label_map, colname='celltype')
workflow.save_seurat()

## Run differential expression externally in R or internally in Python
workflow.get_DE(DE='../Data/DE/DE_hesc_dataset.csv')
```

#### Step 3. Network inference: Either load pre-existing network
```
workflow.get_network(cell_type='embryonic stem cell')
```

#### Or infer network and edge weights
```
## Set up PySCENIC for given dataset
workflow.set_up_pyscenic(species)

## Run PySCENIC (Takes several hours)
workflow.run_pyscenic()

## Post processing of learnt transcriptional network
workflow.learn_weights()
```

#### Step 4. Run PrecICE to identify optimal transcription factor perturbations
```
transition = source_name +'_to_' + target_name
workflow.run_prescient(species='human',
                          network_path=self.network_path,
                          DE_path=workflow.DE_filenames[transition])
```

#### Step 5. Plot results
```
## Plot ranked list of perturbations and relative perturbation magnitude
workflow.perturbation_plot(k=15)

## Plot transcriptional network with applied perturbations
workflow.network_plot()
```
