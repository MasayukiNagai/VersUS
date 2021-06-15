<?php require "templates/header.php"; ?>
<div class="content">

<!-- <h2>Gene Table</h2>

<form method="post">
  <label for="search_id">Search by: </label>
  <select id="search_by" name="search_by">
    <option value="keywords">Keyword</option>
    <option value="gene_symbol">Gene ID</option>
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
    $connection = new PDO($dsn, $username, $password, $options);
    $ec = $_GET['ec'];
    $ec_numbers = explode(".", $ec);
    $current_page = Get('page', 1);
    $num_per_page = 50;
    $start = ($current_page - 1) * $num_per_page;
    if(count($ec_numbers) == 1){
      $condition = "WHERE ec_1 = {$ec_numbers[0]}";
    }
    elseif(count($ec_numbers) == 2){
      $condition = "ec_1 = $ec_numbers[0] 
                    AND ec_2 = {$ec_numbers[1]}";
    }
    elseif(count($ec_numbers) == 3){
      $condition = "ec_1 = $ec_numbers[0] 
                    AND ec_2 = {$ec_numbers[1]}
                    AND ec_3 = {$ec_numbers[2]}";
    }
    elseif(count($ec_numbers) == 4){
      $condition = "ec_1 = {$ec_numbers[0]}
                    AND ec_2 = {$ec_numbers[1]}
                    AND ec_3 = {$ec_numbers[2]}
                    AND ec_4 = {$ec_numbers[3]}";
    }
    else{
      echo "Error: " . $ec_numbers;
    }
    $order = "ORDER BY EC_number";
    $limit = "LIMIT :start, :num_per_page";
    $sql = get_query($condition, $order, $limit);
    // $sql = "SELECT g.gene_id, g.gene_symbol, g.gene_full_name, 
    //                COUNT(m.mutation_id) AS num_vus, 
    //                MAX(m.CADD_score) AS max_cadd, 
    //                g.EC_number
    //         FROM (SELECT * FROM Gene WHERE $condition 
    //               ORDER BY EC_number ASC LIMIT :start, {$num_per_page}) as g
    //         LEFT JOIN Mutation AS m USING(gene_id)
    //         WHERE $condition
    //         GROUP BY g.gene_id";
    $statement = $connection->prepare($sql);
    $statement->bindParam(':start', $start, PDO::PARAM_INT);
    $statement->bindParam(':num_per_page', $num_per_page, PDO::PARAM_INT);
    // $statement->bindParam(':condition', $condition, PDO::PARAM_STR);
    $statement->execute();
    $results = $statement->fetchAll();

    $sql2 = "SELECT COUNT(*) from Gene {$condition}";
    $statement2 = $connection->prepare($sql2);
    $statement2->execute();
    $num_results = $statement2->fetch()[0]; 
    $total_page = ceil($num_results/$num_per_page);
  }catch(PDOException $error) {
    echo $sql . "<br>" . $error->getMessage();
  }
}
else{
    echo "> No results are available.";
}
?>

<?php 
$first_item = ($current_page-1) * $num_per_page + 1;
$last_item = ($current_page * $num_per_page) < $num_results ? $current_page * $num_per_page: $num_results;
?>

<div id="result_header">
  <h2 class="result_count">Items: <?php echo escape($first_item) ?> to <?php echo escape($last_item); ?> of <?php echo escape($num_results) ?> genes</h2>
  <div class="pageforms">
    <form action="#" method="get" class="pageform">
      <button class="page_button" type="submit" name="page" value="1" <?php if($current_page==1){?> disabled="disabled" <?php } ?>>&lt&lt First</button>
      <button class="page_button" type="submit" name="page" value=<?php echo ($current_page-1); ?> <?php if($current_page==1){?>disabled="disabled" <?php }; ?> >&lt Prev</button>
    </form>        
    <form action="#" method="get" class="pageform">
      <span>Page </span><input type="text" id="view_page" name="page" value="<?php echo $current_page; ?>" />/ <span id="total_page"><?php echo $total_page ?></span>
      <input type="submit" id="jump_button" value="Go"/>
    </form>
    <form action="#" method="get" class="pageform"> 
      <button class="page_button" type="submit" name="page" value=<?php echo ($current_page+1); ?> <?php if($current_page==$total_page){?>disabled="disabled" <?php }; ?> >&gt Next</button>     
      <button class="page_button" type="submit" name="page" value=<?php echo $total_page ?> <?php if($current_page==$total_page){?>disabled="disabled" <?php } ?>>&gt&gt Last</button>
    </from>
  </div>
</div>


<?php
$counter = ($current_page-1) * $num_per_page;
if ($results && $statement->rowCount() > 0) { ?>
    <table>
    　<thead>
        <tr>
        <th class="count">#</th> 
        <th class="gene_id">Gene ID</th>
        <th class="enzyme_name">Enzyme Name</th>
        <th class="uniprot_id">Uniprot ID</th>
        <th class="num_vus"># Missense VUS</th>
        <th class="cadd_score">Highest CADD score</th>
        <th class="EC_number">EC #</th>
        </tr>
      </thead>
      <tbody>
      <?php foreach ($results as $row) { 
        $counter += 1; ?>
        <tr>
        <td class="count"><?php echo escape($counter) ?></td>
        <td class="gene_id"><a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>&page=1"><?php echo escape($row["gene_symbol"]); ?></a></td>
        <td class="enzyme_name"><?php echo escape($row["gene_full_name"]); ?></td>
        <td class="uniprot_id"><a href=<?php echo get_uniprot_url($row["uniprot_id"]) ?>><?php echo escape($row["uniprot_id"]) ?></a></td>
        <td class="num_vus"><?php echo escape($row["num_vus"]); ?></td>
        <!-- <td class="cadd_score"><?php echo escape($row["max_cadd"]); ?></td> -->
        <td class="cadd_score"><?php echo escape(number_format((float)$row["max_cadd"], 1, '.', '')); ?></td>
        <td class="EC_number"><?php echo escape($row["EC_number"]); ?></td>
        </tr>
      <?php } ?>
      </tbody>
    </table>
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