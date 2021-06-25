# Get the source data

First, download the data from source for the socio-patterns and the coauthorship:

- [sociopatterns.org/datasets/](http://www.sociopatterns.org/datasets/)
- [github.com/arbenson/ScHoLP-Data](https://github.com/arbenson/ScHoLP-Data)

Place them in the `socio_data` directory

# Format the data

Run the script `format_data.py` to extract distributions and other
representations for each dataset.

# Subsample the dblp coauthorship hypergraph

Run the script `subsample_bfs.py` to have a more reasonable size hypergraph for
simulations.


