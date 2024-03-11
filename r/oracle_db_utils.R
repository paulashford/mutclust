# Oracle db utils
# Ash 27/05/2015

# Ash 20/01/16
# Make code portable by changing classpaths

#install.packages("RJDBC")

connectFGFRlocal <- function(){
  library(RJDBC)
  # NOTE  grant update on FGFR.TBL_ALIGNMENT TO FGFR_READER2;
  #drv <- JDBC("oracle.jdbc.OracleDriver",classPath = "/Users/ash/bioinf/oracle/ojdbc6.jar");
  #drv <- JDBC("oracle.jdbc.OracleDriver",classPath = "/home/ucbtshf/woofgit/fgfr-net-r/shared/ojdbc6.jar");
  drv <- JDBC("oracle.jdbc.OracleDriver",classPath = "../shared/ojdbc6.jar");
  con <- dbConnect(drv, "jdbc:oracle:thin:@localhost:1521:FGFR","xxxxxxxxxxx","xxxxxxx");
  return (con);
 #return(dbListTables(con))
}

connectSINATRA <- function(){
  library(RJDBC);
  library(rJava);
  options(java.parameters="-Xm4g");
  #drv <- JDBC("oracle.jdbc.OracleDriver",classPath = "/Users/ash/bioinf/oracle/ojdbc6.jar");
  drv <- JDBC("oracle.jdbc.OracleDriver",classPath = "../shared/ojdbc6.jar");
  con <- dbConnect(drv, "jdbc:oracle:thin:@sinatra.biochem.ucl.ac.uk:1521:BIOMAPWH","xxxxxxxxxxx","xxxxxx");
  return (con);
}
