# Import main components from their respective modules
from .data_aggregation import aggregate_data
from .visualisation import plot_droplet_analysis

# Import algorithm classes
from .gibbs_dividing_plane import *
from .circle_fitting import *

# Version of the package
__version__ = "0.0.1"

# Define what should be imported with "from contact_angle import *"
__all__ = [
    "aggregate_data",
    "plot_droplet_analysis",
]
