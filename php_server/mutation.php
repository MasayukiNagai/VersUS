<?php require "templates/header.php"; ?>
<div class="content">

<?php 
try{
  require "config.php";
  require "common.php";
  $gene_id = $_GET['gene_id'];
  $current_page_count = $_GET['page'];
  $num_per_page = 50;
  $start = ($current_page_count - 1) * $num_per_page;
  $connection = new PDO($dsn, $username, $password, $options);

  $sql1 = "SELECT gene_name_short, gene_name_full, EC_number
           FROM Gene WHERE gene_id = $gene_id LIMIT 1";
  $statement1 = $connection->prepare($sql1);
  $statement1->execute();
  $value = $statement1->fetch();
  $gene_name_short = $value['gene_name_short'];
  $gene_name_full = $value['gene_name_full'];
  $ec_number = $value['EC_number'];

  $sql2 = "SELECT * FROM Mutation WHERE gene_id = $gene_id 
           ORDER BY CADD_score DESC
           LIMIT :start, $num_per_page";
  $statement2 = $connection->prepare($sql2);
  $statement2->bindParam(':start', $start, PDO::PARAM_INT);
  $statement2->execute();
  $results = $statement2->fetchAll();

  $sql3 = "SELECT COUNT(*) FROM Mutation WHERE gene_id = $gene_id";
  $statement3 = $connection->prepare($sql3);
  $statement3->execute();
  $value3 = $statement3->fetch();
  $total_page_count = ceil($value3[0]/$num_per_page);
  
}catch(PDOException $error) {
  echo $sql2 . "<br>" . $error->getMessage();
}?>

<div id="mutations">

<div id="header_mutation">
  <h1 class="gene_name_short"><?php echo escape($gene_name_short); ?></h1>
  <h2 class="gene_name_full"><?php echo escape($gene_name_full); ?></h2>
</div>


<?php
if ($results && $statement2->rowCount() > 0) { ?>
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
    <?php foreach ($results as $row) { ?>
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

<div class="pagination">
  <?php 
  
  if($current_page_count > 1){?>
  <a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>&page=<?php echo $current_page_count-1 ?>">&laquo;</a>
  <?php }
  $page_count = 1;
  while($page_count <= $total_page_count){?>
  <a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>&page=<?php echo $page_count ?>"><?php echo escape($page_count) ?></a>
  <?php $page_count++;
  } 
  if($current_page_count < $total_page_count){?>
  <a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>&page=<?php echo $current_page_count+1 ?>">&raquo;</a>
  <?php }?> 
</div>

<?php } else { ?>
  > No resultss are available.
<?php } ?>

</div>
</div>
<?php require "templates/footer.php"; ?>