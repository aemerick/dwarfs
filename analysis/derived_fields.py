import yt
import numpy as np
import cgs   as cgs
# Define field functions here

def _CellVolume(field, data):
    return (data['dx'] * data['dy'] * data['dz']).convert_to_units('cm**3')

def _gas_column_density(field,data):
    print " ----- WARNING ----- "
    print "Column calculated using horrible assumptions on\n"
    print "composition and mean molecular weight"   
    print " ----- ------ ------ "

    if data.has_field_parameter('dens'):
        density_name = 'dens'
    else:
        density_name = 'density'

    mu = np.ones(np.size(data['dx']))

    mu[data['Temperature'] < 1.0E4] = 1.21 # assume neutral
    mu[data['Temperature'] > 1.0E4] = 0.61 # assume ionized

    return data['dx'] * data[density_name] / (cgs.mp * mu)


# dictionary of field names
function_dict = {'CellVolume': _CellVolume,
                 'gas_column_density': _gas_column_density}

# display name dictionary
display_name_dict = {'CellVolume': r"\rm{cm}^{3}",
                     'gas_column_density': r'\rm{cm}^{2}'}

# units dictionary
units_dict = {'CellVolume': "cm**(-3)",
              'gas_column_density': "cm**(-2)"}

# field names
fields = function_dict.keys()

def add_field(field):
    """
        Wrapper to yt.add_field function. Supply field name from
        list of possible fields given by derived_fields.fields
    
        Parameters
        ----------
        field : string
            Name of field to add. Must already have function
            definitions contained within the derived_fields.
            See derived_fields.fields for list of field names.
    """

    yt.add_field(field, function = function_dict[field],
                        display_name = display_name_dict[field],
                        units = units_dict[field])

