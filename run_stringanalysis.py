from core.policies.membrane_policy import Membrane
from core.policies.string_policy import String

if __name__=='__main__':
    # Example usage:
    csv_path = 'path/to/your/membrane_data.csv'
    membrane = Membrane(csv_path)
    psii_string = String(membrane)

    # Now you can access psii_strings or other information as needed
    print(psii_string.psii_strings)