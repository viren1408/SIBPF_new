import argparse
import json
from crossmatching_function import Crossmatching #use this for the previous version
from SIBPF import sibpf #new version
# Import the given subroutine here
# Make sure you have the necessary imports and functions defined in the 'Trial_onlyspidx' submodule

# Create a parser
parser = argparse.ArgumentParser(description="Crossmatching and spectral analysis")

# Add an argument for the configuration file
parser.add_argument('--config', required=True, help="Path to configuration JSON file")

# Parse the arguments
args = parser.parse_args()

# Read the configuration from the JSON file
with open(args.config, 'r') as config_file:
    config = json.load(config_file)

# Call the Crossmatching function with the provided arguments from the config
spectral_index, pulsar_candidates = sibpf(
    dir=config['dir'],
    region=config['region'],
    file_path_image=config['image'],
    file_path_spidx=config['spidx'],
    file_path_TGSS=config['tgss'],
    show_matches=config['show_matches'],
    get_spectral_index=config['get_spectral_index'],
    get_candidates=config['get_candidates'],
    get_pulsars=config['get_pulsars']
)

# Print the results if needed
print("Spectral Index:")
print(spectral_index)

print("Pulsar Candidates:")
print(pulsar_candidates)
