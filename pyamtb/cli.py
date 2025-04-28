import argparse
import os
from .tight_binding_model import calculate_band_structure, create_pythtb_model
from .parameters import Parameters
from .read_datas import read_poscar
from .check_distance import calculate_distances

def main():
    parser = argparse.ArgumentParser(description='PyAMTB - Tight-binding model calculations')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands', required=True)
    
    # Model calculation command
    calc_parser = subparsers.add_parser('calculate', help='Calculate tight-binding model')
    calc_parser.add_argument('--config', type=str, help='Path to configuration file')
    calc_parser.add_argument('--poscar', type=str, help='Path to POSCAR file')
    calc_parser.add_argument('--output', type=str, help='Output filename')
    
    # Distance calculation command
    dist_parser = subparsers.add_parser('distance', help='Calculate distances between atoms')
    dist_parser.add_argument('--poscar', type=str, required=True, help='Path to POSCAR file')
    dist_parser.add_argument('--element1', type=str, default="Mn", help='First element type')
    dist_parser.add_argument('--element2', type=str, default="N", help='Second element type')
    
    args = parser.parse_args()
    
    if args.command == 'calculate':
        # Load configuration
        if args.config:
            params = Parameters(args.config)
        else:
            params = Parameters()
            
        # Set output filename if provided
        if args.output:
            params.output_filename = args.output
            
        # Set POSCAR file if provided
        if args.poscar:
            poscar_filename = args.poscar
        else:
            poscar_filename = os.path.join(params.savedir, params.output_filename + ".vasp")
            
        # Create and calculate model
        model = create_pythtb_model(params)
        calculate_band_structure(model, params)
        print(f"Calculation completed! Results saved to {params.output_filename}.{params.output_format}")
        
    elif args.command == 'distance':
        # Calculate distances between atoms
        distances = calculate_distances(args.poscar, args.element1, args.element2)
        print(f"\nFound {len(distances)} distances between {args.element1} and {args.element2} atoms")
        
    else:
        parser.print_help()

if __name__ == '__main__':
    main() 