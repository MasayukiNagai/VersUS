<?php require "templates/header.php"; ?>
<div class="content">

<?php 
try{
  require "config.php";
  require "common.php";
  $gene_id = $_GET['gene_id'];
  $connection = new PDO($dsn, $username, $password, $options);
  $sql1 = "SELECT gene_name_short, gene_name_full, EC_number
           FROM Gene WHERE gene_id = $gene_id LIMIT 1";
  $statement1 = $connection->prepare($sql1);
  $statement1->execute();
  $value = $statement1->fetch();
  $gene_name_short = $value['gene_name_short'];
  $gene_name_full = $value['gene_name_full'];
  $ec_number = $value['EC_number'];
  $sql = "SELECT * FROM Mutation WHERE gene_id = $gene_id LIMIT 50";
  $statement = $connection->prepare($sql);
  $statement->execute();
  $result = $statement->fetchAll();
}catch(PDOException $error) {
  echo $sql . "<br>" . $error->getMessage();
}?>

<div id="mutations">

<div id="header_mutation">
  <h1 class="gene_name_short"><?php echo escape($gene_name_short); ?></h1>
  <h2 class="gene_name_full"><?php echo escape($gene_name_full); ?></h2>
</div>



<?php
if ($result && $statement->rowCount() > 0) { ?>
<table>
  <thead>
    <tr>
      <th class="variation">Missense Variation</th>
      <th class="clinvar_link">Clinvar Link</th>
      <th class="cadd_score">CADD Score</th>
      <th class="gnomAD_AF">gnomAD AF</th>
      <th class="pdb">PDB</th>
    </tr>
  </thead>
  <tbody>
    <?php foreach ($result as $row) { ?>
    <tr>
      <td class="variation"><?php echo escape($row["ref_pos_alt"]); ?></td>
      <td class="clinvar_link"><a href=<?php echo $row["clinvar_link"] ?>>link</a></td>
      <td class="cadd_score"><?php echo escape($row["CADD_score"]); ?></td>
      <td class="gnomAD_AF"><?php echo escape($row["gnomAD_AF"]); ?></td>
      <td class="pdb"><?php echo escape($row["pdb"]); ?></td>
    </tr>
    <?php } ?>
  </tbody>
</table>
  <?php } else { ?>
    > No results are available.
  <?php } ?>

</div>

</div>
<?php require "templates/footer.php"; ?>