
# dpMC
dpMC is developed for prediction of missed cleavage on protein sequences digested by trypsin based on deep-learning.
<div align=center><img src='/pics/start.png'/></div>


## Major functions
* Preprocess fasta file and identified peptides file.
* Train model for missed cleavage prediction for local specific experiments with fine-tuning.
* In silico digestion of fasta file using trained model.


## System requirements & Runtime Environment(RE) Confirmation
* Both of Windows and Linux platforme are supported.
* NVIDIA Graphics Processing Unit (GPU) is highly reconmmended; Central Processing Unit (CPU) calculation is also available but depreciated;
* NVIDIA CUDA 10.0+, CuDNN 7.6.5+ are recomended.
* Keras with Tensorflow backend.

dpMC was developed under Python 3.6.5(Anaconda3 5.2.0 64-bit) with keras tensorflow-gpu backend. Hardware including GPU card NVIDIA GeForce 1080Ti, CPU i7-8086K and 128GB RAM were utilized. 


## Installation 
**1. Installation of Python (Anaconda is recommended)**

   * Anaconda installer can be downloaded from the [Anaconda official site](https://www.anaconda.com/products/individual).
   * Official python installer can be downloaded from the [Python official site](https://www.python.org/downloads/).

**2. Installation of associated packages**

   * Install Tensorflow using `conda install(recommended)` or pip:
   
      * *`conda install -c conda-forge tensorflow`* or *pip install --upgrade tensorflow*
   
   * Install Tensorflow with GPU supported using `conda install(recommended)` or pip:
   
      * *`conda install -c anaconda tensorflow-gpu`* or *pip install --upgrade tensorflow-gpu*
 
   * Install Keras using `conda install(recommended)` or pip:
   
      * *`conda install -c conda-forge keras`* or *pip install keras*

   * Other associated packages including **``os,re,datetime,Bio,pandas,numpy,random,fnmatch``** can also be installed using `conda install(recommended)` or pip.
   

## Files needed for dpMC
* For the training of dpMC model, **peptides files** from two searching software are supported:
   * **Spectronaut**, `the experimental library(.xls) file built by Pulsar in Spectronaut` is supported.
   * **MaxQuant**, `the "peptides.txt" file generated from MaxQuant` is supported.
* **Reference fasta file** for dpMC. For annotating peptides in the training file, corresponding proteome fasta file is needed (for now the standard Uniprot fasta file format is supported).
* **Pretrained model** for fine-tuning. Fine-tuning is provided when training model of dpMC and the pretrained models are provided in the [models](models/) folder.
* **Trained model** for digestion. This file is needed when you have trained your dpMC model and ready to digest the fasta file.
* **Fasta file to be in silico digested**. Customized headers in the fasta file are supported.


## Preprocessing of datasets
For the training of dpMC, only high qualified singly identified missed cleavage sites or cleavae sites are used. Ambiguous cleavage sites which are ovserved having both of cleavaged and missed cleavaged states are abandoned. This function has been integrated in parsing peptides file.

## Procedures to train dpMC model
1) Start using dpMC by opening command interpreter *`cmd.exe`* in windows platform or *`shell`* in Liux platform.
2) Run dpMC by calling *`python`* program: `python dpMC_main.py`.
3) After entering the commond line, follow the prompt and enter `1` to select `training` dpMC.
4) Next, please set your working directory after the prompt which will store all your training models.
5) Put the absolute path of your `reference fasta file` and `pretrained model` after corresponding prompt.
6) Next, please select your identified peptides file either from **MaxQuant** or **Spectronaut**.
7) The trained model can be found under folder `./working directory/dpMC/md/XXX-XX-XX_XX_XX_XX_XXXXXX/`, please keep the best model based on the 'validation loss' for future digestion.

  * The examples for training of dpMC based on **MaxQuant** file:
> <div align=center><img src='/pics/train_MQ.PNG'/></div>


  * The examples for training of dpMC based on **Spectronaut** file:
> <div align=center><img src='/pics/train_SP.PNG'/></div>


## Procedures to digest protein sequences
1) Start using dpMC by opening command interpreter *`cmd.exe`* in windows platform or *`shell`* in Liux platform.
2) Run dpMC by calling *`python`* program: `python dpMC_main.py`.
3) After entering the commond line, follow the prompt and enter `2` to select `digestion` by dpMC.
4) Next, please set your working directory after the prompt which will store all your training models.
5) Put the absolute path of your `fasta file for digestion` and `model for digestion` after corresponding prompt.
6) The result can be found under folder `./working directory/dpMC/Digested/dpMC_fasta_digested.txt`.

  **The progression informaiton will be shown in the progress bar in commond window.**

* The examples for prediction by dpMC:
> <div align=center><img src='/pics/digest.PNG'/></div>


## Notes for the files generated by dpMC

* **The model files**. The model files can be generated during training process by selecting *`train`* function. All the model files are in the following format:
`./dpMC-XXX-Y.YYYYY-Z.ZZZZZ.h5/`, in which `XXX` denotes the *epoch* of the model, `Y.YYYYY` denotes the *training loss* for this epoch of training, `Z.ZZZZZ` denotes the *validation loss* for this epoch of training. 
   * For the `training` function, the model files can be found under directory `./working directory/dpMC/md/XXX-XX-XX_XX_XX_XX_XXXXXX/`.

* **The digestion files**. The digestion files can be generated during prediction of digestions on fasta files by selecting *`predict`* function. The results contain four columns: 1) `ID`: protein ID, 2) `Sequence`: digested peptide sequence, 3) `PEPlength`: peptide length , 4) `Pos_start`: start position in protein sequence for given peptide, 5) `Pos_end`: end position in protein sequence for given peptide, 6) `Pos_KR`: lysine or arginine position in protein sequence for given peptide, 7) `stats_KR`: predicted probabilities of missed cleavages for `Pos_KR`, which the `Pos_KR` less than 8 were assigned as 1 and will not be cleaved, 8) `MC`: the number of missed cleavages.
   * For the `digestion` function, the digested file can be found in `./working directory/dpMC/Digested/ dpMC_fasta_digested.txt`.

## Compiled software of dpMC

**To provide researcher one more handy usage of dpMC, we also provide the compiled software on Windows system, please download [dpMC](https://zenodo.org/record/4592409#.YIbVT2gzaUk) and associated [user manual](https://zenodo.org/record/4592409#.YIbVT2gzaUk) for further information.**

## Contacts
Please submit any issues happened to dpMC on issues section or send Emails to dpMC.sun@gmail.com.
