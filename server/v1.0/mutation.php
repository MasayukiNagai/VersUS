<?php require "templates/header.php"; ?>
<div class="content">

<?php
require "config.php";
require "common.php";
try{
  $connection = new PDO($dsn, $username, $password, $options);
  $gene_id = $_GET['gene_id'];
  $current_page = GetPost('page', 1);
  $num_per_page = 50;
  $start = ($current_page - 1) * $num_per_page;

  $sql1 = "SELECT gene_symbol, gene_full_name, EC_number
           FROM Gene WHERE gene_id = $gene_id LIMIT 1";
  $statement1 = $connection->prepare($sql1);
  $statement1->execute();
  $value = $statement1->fetch();
  $gene_symbol = $value['gene_symbol'];
  $gene_full_name = $value['gene_full_name'];
  $ec_number = $value['EC_number'];

  $sql2 = "SELECT * FROM Mutation WHERE gene_id = $gene_id
           ORDER BY CADD_score DESC
           LIMIT :start, :num_per_page";
  $statement2 = $connection->prepare($sql2);
  $statement2->bindParam(':start', $start, PDO::PARAM_INT);
  $statement2->bindParam(':num_per_page', $num_per_page, PDO::PARAM_INT);
  $statement2->execute();
  $results = $statement2->fetchAll();

  $sql3 = "SELECT COUNT(*) FROM Mutation WHERE gene_id = $gene_id";
  $statement3 = $connection->prepare($sql3);
  $statement3->execute();
  $num_results = $statement3->fetch()[0];
  $total_page = ceil($num_results/$num_per_page);

}catch(PDOException $error) {
  echo $sql2 . "<br>" . $error->getMessage();
}?>

<div id="header_mutation">
  <div id="header_mutation_text">
    <h1 class="gene_symbol"><?php echo escape($gene_symbol); ?></h1>
    <h2 class="gene_full_name"><?php echo escape($gene_full_name); ?></h2>
  </div>
</div>

<?php require "templates/result_header.php" ?>

<?php
$counter = ($current_page-1) * $num_per_page;
if ($results && $statement2->rowCount() > 0) { ?>
<table>
  <thead>
    <tr>
      <th class="count">#</th>
      <th class="variation">Missense Variation</th>
      <th class="clinvar_link">Clinvar Link</th>
      <th class="cadd_score">CADD Score</th>
      <th class="gnomAD_AF">gnomAD AF</th>
      <th class="pdb">PDB</th>
    </tr>
  </thead>
  <tbody>
    <?php foreach ($results as $row) {
      $counter += 1;?>
    <tr>
      <td class="count"><?php echo escape($counter) ?></td>
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
  > No resultss are available.
<?php } ?>

</div>
<?php require "templates/footer.php"; ?>
