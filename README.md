# imits_importer
Repository for the EMBL-EBI Imits importer.

README

The imits restful importer requires a configuration file CFG.pl to be set up locally and excluded from GitHub using gitignore. The contents of the configuration file (shown below for reference) should be set up with the variables shown for local database access and values supplied by Mouse Informatics, EMBL-EBI for $IMITSPASSWD and $IMITSUSER:




package CFG;

# -----------------------------------------------------------------------------
# # # Access to the database
# # # -----------------------------------------------------------------------------
 $HOST       = 'XXXXXXXX';
 $USER       = 'XXXXXXXX';
 $PASSWD     = 'XXXXXXXX';
 $DATABASE   = 'XXXXXXXX';
 $PORT       = 'XXXXXXXX';

 $DSN = "dbi:mysql:$DATABASE:$HOST:$PORT";


 $IMITSPASSWD= 'XXXXXXXX';
 $IMITSUSER = "XXXXXXXX";




Please refer to import_iMits_RESTFUL.pl for usage.

