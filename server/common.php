<?php

/**
  * Escapes HTML for output
  *
  */

function escape($html) {
  return htmlspecialchars($html, ENT_QUOTES | ENT_SUBSTITUTE, "UTF-8");
}

function Get($index, $defaultValue) {
  return isset($_GET[$index]) ? $_GET[$index] : $defaultValue;
}

function GetPost($index, $defaultValue) {
  return isset($_POST[$index]) ? $_POST[$index] : $defaultValue;
}

function get_uniprot_url($uniprot_id) {
  $url = "https://www.uniprot.org/uniprot/";
  return $url . $uniprot_id;
}

function get_alphafold_url($uniprot_id) {
  $url = "https://alphafold.ebi.ac.uk/entry/";
  return $url . $uniprot_id;
}

function get_query($condition, $order, $limit){
  if($condition == ""){
    $sql = "SELECT g.gene_id, g.gene_symbol, g.gene_full_name,
                  g.uniprot_id, g.EC_number,
                  COUNT(m.mutation_id) AS num_vus,
                  MAX(m.CADD_score) AS max_cadd
            FROM Gene AS g
            LEFT JOIN Mutation AS m USING(gene_id)
            GROUP BY g.gene_id
            {$order} {$limit}";
  }
  else{
    $sql = "SELECT g.gene_id, g.gene_symbol, g.gene_full_name,
                  g.uniprot_id, g.EC_number,
                  COUNT(m.mutation_id) AS num_vus,
                  MAX(m.CADD_score) AS max_cadd
            FROM (SELECT * FROM Gene {$condition}) AS g
            LEFT JOIN Mutation AS m USING(gene_id)
            GROUP BY g.gene_id
            {$order} {$limit}";
  }
  return $sql;
}
