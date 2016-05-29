## ModelWork - MATLAB suite for model building and model fitting

This suite contains a set of standardized functions for model building and model evaluation (model fitting via maximum likelihood and MCMC sampling, and computation of model comparison metrics).

I started creating this during my PhD, and it is mostly for personal usage.
One day I might clean up the code enough to be actually usable by other people.

The suite name is from the sentence *"please model work"* that all modellers utter before running a new model, or to an old model after a long debugging session. This etymology emerged during a chat I had with [Will Adler](https://github.com/wtadler) and it is now canon.

---

### Setting up the pipeline for a new project

#### Convert the raw datasets to standard format

Write a function `project_convertData.m` that converts the original raw data files into the standardized **ModelWork** format. This consists of a single `data` cell array, with one element per subject. Each element has fields:
- `id`: the subject's identificative number.
- `name`: the subject's identificative string.
- `X`: a struct that contains the standardized data. `X` can have several subfields, e.g. divided by task type. The data matrices are usually contained in cell arrays, with one data matrix per noise condition (where applicable). The first column of each data matrix is usually reserved for the trial number or session number.
- Additional fields might contain other useful information.

#### Write the log likelihood function

Write a function `project_like.m` that computes the log likelihood of a dataset. The function takes as input:
- A struct `X` of data as defined above.
- A model-parameter structure `mp` (see below). 
- Auxiliary struct `infostruct`. This struct usually contains bulky auxiliary data used during computations (e.g., a lookup table); `infostruct` is precomputed or loaded from file at the beginning of the model fitting procedure and then passed around.
- A `flags` array (rarely used).

The log likelihood function will usually loop over conditions and calculate the log likelihood for each condition by calling other condition- and model- dependent subroutines. Persistent variables should be used to prevent recomputation of conditions whose parameters have not changed since the last call to *project_like.m*.

#### Write the model and parameter setup function

Write a function `project_setupModel.m` that sets up the various models and their parameters. Use one of the existing files as template.

#### Write auxiliary functions

Knowing the structure of the models and parameters, you can now write a bunch of auxiliary functions:

- `project_logPrior.m` computes the log prior for each parameter (depending on the model).
- `project_getModelName.m` returns the model string associated with a given model vector.
