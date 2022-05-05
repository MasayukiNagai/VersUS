<?php

/**
  * Configuration for database connection
  */

$server     = "server";
$username   = "username";
$password   = "password";
$db_name    = "db_name";
$dsn        = "mysql:host=$server;dbname=$db_name;port=3306;";
$options    = array(
                PDO::ATTR_ERRMODE => PDO::ERRMODE_EXCEPTION
              );
