# Declaring global variables for dplyr and data.table column names
utils::globalVariables(c("sampleid","treatmentid",'read.csv','samplename','.',
                         'Symbol', 'feature', 'Control', 'Low', 'Middle',
                         'High', 'verbose', 'dose_level', 'individual_id',
                         'duration_h', 'viability', '.SD', 'durations',
                         'tSetName', '.intern', 'controlLevels', 'treatmentLevels',
                         'dose', 'compound', 'duration'))