<?php require "templates/header.php"; ?>
<h2>Mutation Table</h2>

<?php 
try{
  require "config.php";
  require "common.php";
  $ec_number = $_GET['ec_number'];
  $connection = new PDO($dsn, $username, $password, $options);
  $sql1 = "SELECT gene_id, gene_name_short FROM Gene WHERE EC_number = :ec_number LIMIT 1";
  $statement1 = $connection->prepare($sql1);
  $statement1->bindParam(':ec_number', $ec_number, PDO::PARAM_STR);
  $statement1->execute();
  $value = $statement1->fetch();
  $gene_id = $value['gene_id'];
  $gene_name = $value['gene_name_short'];
  ?>
<h3><?php echo escape($gene_name); ?></h3>
<?php
  $sql = "SELECT * FROM Mutation WHERE gene_id = $gene_id LIMIT 50";
  $statement = $connection->prepare($sql);
  $statement->execute();
  $result = $statement->fetchAll();
}catch(PDOException $error) {
  echo $sql . "<br>" . $error->getMessage();
}
?>

<?php
if ($result && $statement->rowCount() > 0) { ?>
    <table>
      <thead>
<tr>
  <th>#</th>
  <th>Missense Variation</th>
  <th>Clinvar Link</th>
  <th>CADD Score</th>
  <th>gnomAD AF</th>
  <th>PDB</th>
</tr>
      </thead>
      <tbody>
  <?php foreach ($result as $row) { ?>
      <tr>
<td><?php echo escape($row["mutation_id"]); ?></td>
<td><?php echo escape($row["ref_pos_alt"]); ?></td>
<td><a href=<?php echo $row["clinvar_link"] ?>>link</a></td>
<td><?php echo escape($row["CADD_score"]); ?></td>
<td><?php echo escape($row["gnomAD_AF"]); ?></td>
<td><?php echo escape($row["pdb"]); ?></td>
      </tr>
    <?php } ?>
      </tbody>
  </table>
  <?php } else { ?>
    > No results are available.
  <?php } ?>