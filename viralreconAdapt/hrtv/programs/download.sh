mkdir input
mkdir water_input

project="HRTV"

dx select "$project"

dir="2023-01-31"

dx download "$dir"/1_raw_data/*.gz -o ./input/

mv ./input/Water* ./water_input/