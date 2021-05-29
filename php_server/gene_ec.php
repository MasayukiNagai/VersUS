<?php require "templates/header.php"; ?>
<div class="content">

<!-- <h2>Gene Table</h2>

<form method="post">
  <label for="search_id">Search by: </label>
  <select id="search_by" name="search_by">
    <option value="keywords">Keyword</option>
    <option value="gene_name_short">Gene ID</option>
    <option value="uniprot_id">Uniprot ID</option>
  </select>
  <input type="text" id="keyword" name="keyword">
  <input type="submit" name="submit" value="Search">
</form> -->

<?php 
if(!isset($_POST['submit'])){
  try{
    require "config.php";
    require "common.php";
    $ec = $_GET['ec'];
    $ec_numbers = explode(".", $ec);
    if(count($ec_numbers) == 1){
      $condition = "ec_1 = $ec_numbers[0]";
    }
    elseif(count($ec_numbers) == 2){
      $condition = "ec_1 = $ec_numbers[0] 
                    AND ec_2 = $ec_numbers[1]";
    }
    elseif(count($ec_numbers) == 3){
      $condition = "ec_1 = $ec_numbers[0] 
                    AND ec_2 = $ec_numbers[1]
                    AND ec_3 = $ec_numbers[2]";
    }
    elseif(count($ec_numbers) == 4){
      $condition = "ec_1 = $ec_numbers[0] 
                    AND ec_2 = $ec_numbers[1]
                    AND ec_3 = $ec_numbers[2]
                    AND ec_4 = $ec_numbers[3]";
    }
    else{
      echo "Error: " . $ec_numbers;
    }
    $current_page_count = $_GET['page'];
    $num_per_page = 50;
    $start = ($current_page_count - 1) * $num_per_page;
    $connection = new PDO($dsn, $username, $password, $options);
    $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                   COUNT(m.mutation_id) AS num_vus, 
                   MAX(m.CADD_score) AS max_cadd, 
                   g.EC_number
            FROM Gene as g
            LEFT JOIN Mutation AS m USING(gene_id)
            WHERE $condition
            GROUP BY g.gene_id
            ORDER BY g.EC_number ASC
            LIMIT :start, $num_per_page";
    $statement = $connection->prepare($sql);
    $statement->bindParam(':start', $start, PDO::PARAM_INT);
    // $statement->bindParam(':condition', $condition, PDO::PARAM_STR);
    $statement->execute();
    $result = $statement->fetchAll();

    $sql2 = "SELECT COUNT(*) from Gene WHERE $condition";
    $statement2 = $connection->prepare($sql2);
    $statement2->execute();
    $result2 = $statement2->fetch();
    $total_page_count = ceil($result2[0]/$num_per_page);
  }catch(PDOException $error) {
    echo $sql . "<br>" . $error->getMessage();
  }
}
else{
    echo "> No results are available.";
}
?>

<?php
if ($result && $statement->rowCount() > 0) { ?>
<table>
　<thead>
    <tr>
      <!-- <th>#</th> -->
      <th class="gene_id">Gene ID</th>
      <th class="enzyme_name">Enzyme Name</th>
      <th class="num_vus"># Missense VUS</th>
      <th class="cadd_score">Highest CADD score</th>
      <th class="EC_number">EC #</th>
    </tr>
  </thead>
  <tbody>
  <?php foreach ($result as $row) { ?>
    <tr>
      <!-- <?php print_r(array_keys($row))?> -->
      <td class="gene_id"><a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>"><?php echo escape($row["gene_name_short"]); ?></a></td>
      <td class="enzyme_name"><?php echo escape($row["gene_name_full"]); ?></td>
      <td class="num_vus"><?php echo escape($row["num_vus"]); ?></td>
      <td class="cadd_score"><?php echo escape($row["max_cadd"]); ?></td>
      <td class="EC_number"><?php echo escape($row["EC_number"]); ?></td>
    </tr>
    <?php } ?>
      </tbody>
</table>
  
<!-- <p><?php print_r(array_values($result2));?></p>
<p><?php echo escape($total_page_count);?></p> -->


<div class="pagination">
  <?php 
  if($current_page_count > 1){?>
  <a href="gene_ec.php?ec=<?php echo $ec ?>&page=<?php echo $current_page_count-1 ?>">&laquo;</a>
  <?php }
  $page_count = 1;
  while($page_count <= $total_page_count){?>
  <a href="gene_ec.php?ec=<?php echo $ec ?>&page=<?php echo $page_count ?>"><?php echo escape($page_count) ?></a>
  <?php $page_count++;
  } 
  if($current_page_count < $total_page_count){?>
  <a href="gene_ec.php?ec=<?php echo $ec ?>&page=<?php echo $current_page_count+1 ?>">&raquo;</a>
  <?php }?> 
</div>

  <?php } else { ?>
    <p>> No results are available.</p>
  <?php } ?>

<?php
if (isset($_POST['submit'])){?>
  <a href="gene.php">Reset</a>
<?php } ?>

</div>
<?php require "templates/footer.php"; ?>