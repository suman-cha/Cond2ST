# Cond2ST: General Frameworks for Conditional Two-Sample Testing

This repository contains the R codes for reproducing the experimental results presented in the paper "General Frameworks for Conditional Two-sample Testing" [link]. 
The paper introduces general methodologies for conditional two-sample testing, aiming to determine whether two populations have the same distributions after accounting for confounding variables. We introduce two approaches: (1) converting conditional independence tests to conditional two-sample tests and (2) density ratio-based testing. If you are planning to do tests for your own projects, I recommend to use R package named "Cond2ST". 

## Getting Started 

### Installation

To get started with the repository, clone it to your local machine:

```sh
git clone https://github.com/suman-cha/Cond2ST.git
cd Cond2ST
```

### Run Experiments

The experiments in the paper can be reproduced using the R codes provided in the experiments/ directory. Each code corresponds to an experiment discussed in the paper.

1. **Prepare the Data**: The superconductivity dataset used in the real data experiments are included in the `real_examples/data/` folder. Make sure the data is in the correct format before running the scripts. Other simulation data are generated in each code. 

2. **Run Scripts**: Navigate to the `experiments/` folder and run the desired experiment script. For example:
   ```sh
   Rscript experiments/scenarios/scenario1u.R
   ```
   The results will be saved in the `results/` folder.



## Citation

If you find this code useful for your research, please cite the paper:

TBD

## Contact

For questions or issues regarding the code, please do not hesitate to contact us!

Seongchan Lee* : statchan1106@yonsei.ac.kr

Suman Cha* : oldrain123@yonsei.ac.kr 

Ilmun Kim : ilmun@yonsei.ac.kr 

## License

TBD

## Acknowledgments

TBD
