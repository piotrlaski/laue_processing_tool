# Processing Script for Time-Resolved X-Ray Diffraction Data Analysis

An automated shell script for processing time-resolved X-ray diffraction data from pump-probe experiments, handling the complete workflow from raw data export to reflection intensity analysis.

## Getting Started

Clone the repository into your local space and ensure you have the required dependencies installed. Configure the paths and compound settings in the script header, then run `processing.sh` with the appropriate processing option.

## Required Input Files

You will need the following files in a subdirectory named after your compound:

1. **laser_01.inp** - Configuration file for laser-on measurements
2. **dark_01.inp** - Configuration file for laser-off (dark) measurements
3. **[COMPOUND]_mono.h5** - Monochromated diffraction data file
4. **[COMPOUND].mda** - Beam parameter data file

## Usage

Run the script with different processing options:

```bash
./processing.sh [option]
```

### Available Processing Options:

- **c** - Creating directories and .inp files
- **e** - Exporting data
- **s** - Signal processing
- **i** - Integrating reflections
- **d** - Dumping to LaueUtil files AND deleting exported files
- **a** - Complete workflow: exporting, signal processing, integrating, and dumping to LU
- **p** - Plotting integrated powerscan data
- **m** - Matching pairs in dark datasets
- **o** - Orientation matrix refinement
- **hkl** - HKL assignment and plotting Euler angles
- **r** - HKL index to ratio conversion
- **sortav** - Sortav processing for reflection averaging (NOTE: this function is implemented better in the TRL_Mapper project, which merges data and creates PD maps)
- **b** - Save masks
- **u** - Model refinement
- **h** - Display help information

## Configuration

Configure the following paths in the script header:

```bash
COMPOUND="Your_Compound_Name"
PATH_DATA="path/to/your/data"
PATH_SCRIPTS="path/to/scripts"
XTAL_DIR="path/to/xtal"
LAUE_DIR="path/to/laueutil"
```

## Dependencies

- **LaueUtil** - For Laue diffraction data processing
- **XTAL**/**LaueProc** - GPU-accelerated photocrystallographic analysis tools
- Standard Unix tools (bash, awk, bc)

LaueUtil most relevant paper:
Kalinowski, J. A., Fournier, B., Makal, A. & Coppens, P. (2012). J. Synchrotron Rad. 19, 637-646
https://doi.org/10.1107/S0909049512022637

LaueProc most relevant paper:
Szarejko, D., Kamiński, R., Łaski, P. & Jarzembska, K. N. (2020). J. Synchrotron Rad. 27, 405-413
https://doi.org/10.1107/S1600577520000077

## Output

The script processes data through various stages:

1. **Data Export** - Extracts diffraction patterns from raw data
2. **Signal Processing** - Applies filters and corrections
3. **Integration** - Calculates reflection intensities
4. **Pair Matching** - Matches laser-on/laser-off reflection pairs for ratio calculation
5. **Ratio Analysis** - Generates intensity ratios for temperature analysis

Results are organized in structured output directories with processed data ready for further analysis.

## File Structure

- `processing.sh` - Main processing script
- `merger2_02.sh` - Script for merging reflection data from multiple measurements
- `templates/` - Template files for various processing steps

## Advanced Features

- **Automatic pair matching** with adjustable intensity thresholds for dark measurements
- **Batch processing** for multiple datasets
- **Integration with crystallographic refinement tools**
- **Automated quality control** and data validation

For any additional help or bug reports, please contact your research group or submit an issue to the repository.
