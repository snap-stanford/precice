# PreciCE

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

#### Step 1. Initialize PreciCE data processing workflow and preprocess dataset
```
path = '../data/Friedman.h5ad'
raw_adata = sc.read_h5ad('../data/Friedman.h5ad')
workflow = precice(adata=raw_adata, 
                   path = path, 
                   cell_filter=True)
```

#### Step 2. Run batch correction (optional)
```
workflow.set_up_scvi(batch_key='day')

## Visualize effects

workflow.scvi_plot_setup()
sc.pl.umap(workflow.adata, color='day')
```

#### Step 3. Compute differential expression
```
### Identify source and target cells
source_name = 'stem'
target_name = 'meso'

workflow.get_DE(source_name=source_name, target_name=target_name)
```

#### Step 4. Network inference: Either load pre-existing network
```
workflow.get_network(cell_type='embryonic stem cell')
```

#### Or infer network and edge weights [Can take several hours]
```
## Set up PySCENIC for given dataset
workflow.set_up_pyscenic(species)

## Run PySCENIC (Takes several hours)
workflow.run_pyscenic()

## Post processing of learnt transcriptional network
workflow.learn_weights()
```

#### Step 5. Run PrecICE to identify optimal transcription factor perturbations
```
transition = source_name +'_to_' + target_name
python_path = '/user/bin/python'
workflow.run_precice(species='human', python_path=python_path,
                     network_path=workflow.network_path,
                     DE_path=workflow.DE_filenames[transition])
```

#### Step 6. Plot results
```
## Plot ranked list of perturbations and relative perturbation magnitude
workflow.perturbation_plot(k=15)

## Plot transcriptional network with applied perturbations
workflow.network_plot()
```
