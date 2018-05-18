When used for the pipeline, this directory is expected to contain at least the following files:
* `jdkLogConfig.properties`: a JAVA property file allowing to configure logging messages using the package `java.util.logging`
* `log4j.properties`: a JAVA property file allowing to configure logging messages using the library log4j
* `log4j2.xml`: a XML file allowing to configure logging messages using the library log4j2
* a fat JAR of the Bgee pipeline, shipping all its dependencies

The default configuration of `java.util.logging` will be overridden using `jdkLogConfig.properties`.  
If you do not want this behavior to take place, use the system property `java.util.logging.config.file` to point to a different configuration file.
  
The default logging level of the SLF4J SimpleLogger will be set to `error`. Use the system property `org.slf4j.simpleLogger.defaultLogLevel` to change this value.
  
The configuration of log4j and log4j2 are simply retrieved from the associated configuration files.

This will cover all loggers used by the dependencies of Bgee (SLF4J SimpleLogger, JDK java.util.logging package, log4, log4j2).

### System properties specific to Bgee
* `bgee.dao.jdbc.username`: username to connect to the Bgee database using JDBC
* `bgee.dao.jdbc.password`: password to connect to the Bgee database using JDBC
* `bgee.dao.jdbc.driver.names`: comma-separated list of classnames of the JDBC Driver implementations to use, e.g.: `com.mysql.jdbc.Driver,net.sf.log4jdbc.sql.jdbcapi.DriverSpy`
* `bgee.dao.jdbc.url`: the JDBC connection URL; it should NOT contain the username and password
