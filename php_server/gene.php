<?php require "templates/header.php"; ?>
<div class="content">

<!-- <h2>Gene Table</h2> -->

<!-- <form method="post">
  <label for="search_id">Search by: </label>
  <select id="search_by" name="search_by">
    <option value="gene_name_short">Gene ID</option>
    <option value="uniprot_id">Uniprot ID</option>
    <option value="keywords">Keyword</option>
  </select>
  <input type="text" id="keyword" name="keyword">
  <input type="submit" name="submit" value="Search">
</form> -->

<?php 
if(!isset($_POST['submit'])){
  try{
    require "config.php";
    require "common.php";
    $current_page_count = $_GET['page'];
    $num_per_page = 50;
    $start = ($current_page_count - 1) * $num_per_page;
    $connection = new PDO($dsn, $username, $password, $options);
    $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                   COUNT(m.mutation_id) AS num_vus, 
                   MAX(m.CADD_score) AS max_cadd, 
                   g.EC_number
            FROM (SELECT * FROM Gene 
                  ORDER BY gene_name_short ASC LIMIT :start, $num_per_page) as g
            LEFT JOIN Mutation AS m USING(gene_id)
            GROUP BY g.gene_id";
    $statement = $connection->prepare($sql);
    $statement->bindParam(':start', $start, PDO::PARAM_INT);
    $statement->execute();
    $results = $statement->fetchAll();

    $sql2 = "SELECT COUNT(*) FROM Gene";
    $statement2 = $connection->prepare($sql2);
    $statement2->execute();
    $result2 = $statement2->fetch();
    $total_page_count = ceil($result2[0]/$num_per_page);
  }catch(PDOException $error) {
    echo $sql . "<br>" . $error->getMessage();
  }
}
else{
  try{
    require "config.php";
    require "common.php";
    $current_page_count = $_GET['page'];
    $num_per_page = 50;
    $start = ($current_page_count - 1) * $num_per_page;
    $connection = new PDO($dsn, $username, $password, $options);
    $search_by = $_POST['search_by'];
    $keyword = $_POST['keyword'];
    if ($search_by == 'gene_name_short'){
      // search by Gene ID
      $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                     COUNT(m.mutation_id) AS num_vus, 
                     MAX(m.CADD_score) AS max_cadd, 
                     g.EC_number
              FROM (SELECT * FROM Gene WHERE gene_name_short LIKE :keyword) AS g
              LEFT JOIN Mutation AS m USING(gene_id)
              GROUP BY g.gene_id
              LIMIT 50";
      $keyword = "$keyword%";
      $sql2 = "SELECT COUNT(*) FROM Gene WHERE gene_name_short LIKE :keyword";
    }
    elseif($search_by == 'uniprot_id'){
      // search by Uniprot ID
      $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                     COUNT(m.mutation_id) AS num_vus, 
                     MAX(m.CADD_score) AS max_cadd, 
                     g.EC_number
              FROM (SELECT * FROM Gene WHERE uniprot_id = :keyword) AS g
              LEFT JOIN Mutation AS m USING(gene_id)
              GROUP BY g.gene_id
              LIMIT 50";
      $sql2 = "SELECT COUNT(*) FROM Gene WHERE uniprot_id = :keyword";
    }
    else{ // search by keywords
      $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                     COUNT(m.mutation_id) AS num_vus, 
                     MAX(m.CADD_score) AS max_cadd, 
                     g.EC_number
              FROM (SELECT * FROM Gene
                    WHERE gene_name_short LIKE :keyword
                       OR gene_name_full LIKE :keyword
                       OR EC_number LIKE :keyword) AS g
              LEFT JOIN Mutation AS m USING(gene_id)
              GROUP BY g.gene_id
              LIMIT 50";
      $keyword = "%$keyword%";
      $sql2 = "SELECT COUNT(*) FROM Gene WHERE gene_name_short LIKE :keyword
                                         OR gene_name_full LIKE :keyword
                                         OR EC_number LIKE :keyword";
    }
    $statement = $connection->prepare($sql);
    $statement->bindParam(':keyword', $keyword, PDO::PARAM_STR);
    $statement->execute();
    $results = $statement->fetchAll();
    
    $statement2 = $connection->prepare($sql2);
    $statement2->bindParam(':keyword', $keyword, PDO::PARAM_STR);
    $statement2->execute();
    $value2 = $statement2->fetch();
    $total_page_count = ceil($value2[0]/$num_per_page);
  }catch(PDOException $error) {
    echo $sql . "<br>" . $error->getMessage();
  }
}

?>

<?php
if ($results && $statement->rowCount() > 0) { ?>
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
  <?php foreach ($results as $row) { ?>
    <tr>
      <!-- <?php print_r(array_keys($row))?> -->
      <td class="gene_id"><a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>&page=1"><?php echo escape($row["gene_name_short"]); ?></a></td>
      <td class="enzyme_name"><?php echo escape($row["gene_name_full"]); ?></td>
      <td class="num_vus"><?php echo escape($row["num_vus"]); ?></td>
      <td class="cadd_score"><?php echo escape($row["max_cadd"]); ?></td>
      <td class="EC_number"><?php echo escape($row["EC_number"]); ?></td>
    </tr>
    <?php } ?>
      </tbody>
</table>
  

<div class="pagination">
  <?php 
  
  if($current_page_count > 1){?>
  <a href="gene.php?page=<?php echo $current_page_count-1 ?>">&laquo;</a>
  <?php }
  $page_count = 1;
  while($page_count <= $total_page_count){?>
  <a href="gene.php?page=<?php echo $page_count ?>"><?php echo escape($page_count) ?></a>
  <?php $page_count++;
  } 
  if($current_page_count < $total_page_count){?>
  <a href="gene.php?page=<?php echo $current_page_count+1 ?>">&raquo;</a>
  <?php }?> 
</div>

<?php } else { ?>
  <p>> No results are available.</p>
<?php } ?>

<?php
if (isset($_POST['submit'])){?>
  <a href="gene.php?page=1">Reset</a>
<?php } ?>

</div>
<?php require "templates/footer.php"; ?>