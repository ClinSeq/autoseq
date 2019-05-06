from reportgen.reporting.util import create_sql_session
from reportgen.reporting.metadata import query_database
from referralmanager.cli.models.referrals import AlasccaBloodReferral, AlasccaTissueReferral


def get_hospital_code(sample_barcode, referral_type, db_config_file):
    """
    Extract the hospital code for the given sample barcode from the referral database

    :param sample_barcode: The specific sample identifier
    :param referral_type: The class that holds info on the columns for the specific referral type in the database table
    :param db_config_file: File holding configurations and credentials for the referral database
    :return: An integer hospital code
    """
    
    # Retrieve the referral information for the specified sample from the database
    session = create_sql_session(db_config_file)
    ref = query_database(sample_barcode, referral_type, session)
    
    # Check that an integer hospital code is found
    if not isinstance(ref.hospital_code, int):
        raise TypeError("Hospital code of type {} found for sample barcode {}, referral type {} and config file {}. Expected integer.".
            format(type(ref.hospital_code), sample_barcode, referral_type.__class__.__name__, db_config_file))

    return ref.hospital_code


def get_hospital_country(hospital_code, hospitalhundreds2country):
    """
    Get the country where the hospital is located based on the hospital code
    
    :param hospital_code: Hospital code from referral
    :param hospitalhundreds2country: Dictionary over which hospital code hundred serie that code for which country
    :return: A string value for the country
    """
    
    if not isinstance(hospital_code, int):
        raise TypeError("get_country() arg 1 (hospital_code) must be an integer. {} supplied.".
                        format(type(hospital_code)))
    
    # Round down to the closest hundred series
    hundred_serie = hospital_code//100*100
    
    # Find the country
    try:
        country = hospitalhundreds2country[hundred_serie]
    except KeyError:
        raise KeyError("No known country for hospital code {}, hundred series {}. Known combinations are: {}".
                       format(hospital_code, hundred_serie, hospitalhundreds2country))
    
    return country


def get_report_type(country, country2reporttype):
    """
    Figure out which report type that should be used for the final pdf, depending on the origin country of a sample

    :param country: The origin country of a sample
    :param country2reporttype: Dictionary over which kind of report each country should have
    :return: A string value describing the report type
    """
    
    # Get the report type
    try:
        report_type = country2reporttype[country]
    except KeyError, e:
        raise KeyError("No known report type for country {}. Known combinations are: {}".
                       format(country, country2reporttype))
    
    return report_type


def only_alascca_class_report(blood_barcode, tissue_barcode, db_config_file):
    """
    Figure out whether only the alascca class results or all the results should be printed on the pdf report,
    depending on which country the sample comes from.

    :param blood_barcode: The specific sample identifier for the blood sample
    :param tissue_barcode: The specific sample identifier for the tissue sample
    :param db_config_file: File holding configurations and credentials for the referral database
    :return: True if only the alascca class is to be reported, False otherwise
    """

    # FIXME: Hard coding these dictionaries here for now. Could ultimately be provided on command line for more flexibility
    # Dictionary over which hospital code hundred series that code for which country
    hospitalhundreds2country = {300: "SWEDEN",
                                400: "NORWAY",
                                500: "DENMARK",
                                600: "FINLAND"}
    
    # Dictionary over which kind of report each country should have
    country2reporttype ={"DENMARK": "alascca_only",
                         "NORWAY": "full",
                         "SWEDEN": "full",
                         "FINLAND": "full"}
    
    # Get the hospital codes from the referrals
    hospital_blood = get_hospital_code(blood_barcode, AlasccaBloodReferral, db_config_file)
    hospital_tissue = get_hospital_code(tissue_barcode, AlasccaTissueReferral, db_config_file)

    # Get the countries from the hospital codes
    country_blood = get_hospital_country(hospital_blood, hospitalhundreds2country)
    country_tissue = get_hospital_country(hospital_tissue, hospitalhundreds2country)
    
    # Get the report types from the countries
    report_type_blood = get_report_type(country_blood, country2reporttype)
    report_type_tissue = get_report_type(country_tissue, country2reporttype)
    
    # Return true (i.e. only the alascca class should be on the pdf report) if both blood and tissue
    # have type "alascca_only"
    if report_type_blood == "alascca_only" and report_type_tissue == "alascca_only":
        return True
    else:  # Otherwise returning false, i.e. alscca_only report should not be used
        return False
