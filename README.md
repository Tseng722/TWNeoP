# TWNeoP: A Framework for Neoantigen Prioritization

TWNeoP is a LightGBM-based framework for prioritizing immunogenic neoantigens in personalized cancer vaccines and infectious disease immunotherapies. It integrates features like binding affinity, hydrophobicity, and sequence homology to predict peptide–MHC (pMHC) immunogenicity. Enhanced by ensemble learning with DeepHLApan and BigMHC scores, TWNeoP is trained on curated datasets (IEDB, PRIME, DeepNeo), achieving an AUC of 0.88 and a positive predictive value of ~0.9 on an independent test set. Case studies on tumor-associated neoantigens and tuberculosis vaccine candidates demonstrate its versatility. Open-source and scalable, TWNeoP streamlines vaccine development and accelerates precision immunotherapy.

## Dependencies
Please make sure the following tools have been installed in your envs.
<pre><code>
conda install -y make gcc_linux-64 gxx_linux-64
git clone https://github.com/karchinlab/bigmhc.git
</code></pre>
python 3.9
pip



## Installation
<pre><code>pip -r requirements.txt</code></pre>
**Or**
<pre><code>conda env create -f environment.yml</code></pre>

## Command Line Usage 
After installation, run the CLI tool:
Parameters
<pre><code>
-i : input txt file
-o : output dictionary
-n : sample name
</code></pre>

<pre><code>
twneop -i <input_file> -o <output_dir> -n <run_name>
</code></pre>
### Example
<pre><code>
python [PATH]/TSA_immuno_prediction/tsa_immuno_prediction/libs/main.py \
-i [PATH]/myfilecopy3.txt \
-o [PATH]/work1791/cindy2270/web/TSA_immuno_prediction/test/test_file/ \
-n 02
</code></pre>

## Author
Developed by Yu Hsuan Tseng as part of tumor neoantigen research pipeline.