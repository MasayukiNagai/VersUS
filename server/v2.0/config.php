<?php

/**
  * Configuration for database connection
  */

$server     = "alucard.csc.depauw.edu";
$username   = "test_user";
$password   = "test_password";
$db_name    = "versus_db2";
$dsn        = "mysql:host=$server;dbname=$db_name;port=3306;";
$options    = array(
                PDO::ATTR_ERRMODE => PDO::ERRMODE_EXCEPTION
              );
