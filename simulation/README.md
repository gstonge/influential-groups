# Hypergraph data

## Pickle file

The hypergraph data are already available in the `socio_data` directory.
They are already under the right format to run simulations.
They are binary files that can be imported using pickle.

Let us know if you encounter any issue.

## Get the source data

Otherwise, you can get the data from source and format them.

First, download the data from source for the socio-patterns and the coauthorship:

- [sociopatterns.org/datasets/](http://www.sociopatterns.org/datasets/)
- [github.com/arbenson/ScHoLP-Data](https://github.com/arbenson/ScHoLP-Data)

The data must be a list of groups.

Place them in the `socio_data` directory

### Format the data

Run the script `format_data.py` to extract distributions and other
representations for each dataset.

### Subsample the dblp coauthorship hypergraph

Run the script `subsample_bfs.py` to have a more reasonable size hypergraph for
simulations.


# Simulation and AME solutions

To run simulation, you need to define the experience parameters in a
configuration file.

Configuration files for the experiences we did are already present in the `dat`
directory.

For instance, to run Monte-Carlo simulation on unix system:

`python simulation.py --dir dat/ --conf Lyon_exp1_sim_orig`

This will run the experience related to configuration file
`Lyon_exp1_sim_orig.json` in the `dat` directory and it will save the results
in the same directory using pickle.

Similarly, to obtain the ame solutions:

`python ame_sol.py --dir dat/ --conf Lyon_exp1_ame`


