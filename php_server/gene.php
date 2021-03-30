<?php require "templates/header.php"; ?>
<h2>Gene Table</h2>

<form method="post">
  <label for="search_id">Search by: </label>
  <select id="search_by" name="search_by">
    <option value="gene_name_short">Gene ID</option>
    <option value="uniprot_id">Uniprot ID</option>
    <option value="keywords">Keywords</option>
  </select>
  <input type="text" id="keyword" name="keyword">
  <input type="submit" name="submit" value="Search">
</form>

<?php 
if(!isset($_POST['submit'])){
  try{
    require "config.php";
    require "common.php";
    $connection = new PDO($dsn, $username, $password, $options);
    $sql = "SELECT * FROM Gene LIMIT 50";
    $statement = $connection->prepare($sql);
    $statement->execute();
    $result = $statement->fetchAll();
  }catch(PDOException $error) {
    echo $sql . "<br>" . $error->getMessage();
  }
}
else{
  try{
    require "config.php";
    require "common.php";
    $connection = new PDO($dsn, $username, $password, $options);
    $search_by = $_POST['search_by'];
    $keyword = $_POST['keyword'];
    if ($search_by == 'gene_name_short'){
      // search by Gene ID
      $sql = "SELECT * FROM Gene WHERE $search_by = :keyword LIMIT 50";
    }
    elseif($search_by == 'uniprot_id'){
      // search by Uniprot ID
      $sql = "SELECT * FROM Gene WHERE $search_by = :keyword LIMIT 50";
    }
    else{ // search by keywords
      $sql = "SELECT * FROM Gene WHERE (gene_name_short LIKE :keyword
                                         OR gene_name_full LIKE :keyword)
                                        LIMIT 50";
      $keyword = "%$keyword%";
    }
    $statement = $connection->prepare($sql);
    // $statement->bindColumn(':search_by', $search_by, PDO::PARAM_STR);
    $statement->bindParam(':keyword', $keyword, PDO::PARAM_STR);
    $statement->execute();
    $result = $statement->fetchAll();
    // echo $sql . "<br>" .$search_by .$keyword. "<br>" ;
  }catch(PDOException $error) {
    echo $sql . "<br>" . $error->getMessage();
  }
}

?>

<?php
if ($result && $statement->rowCount() > 0) { ?>
    <table>
      <thead>
<tr>
  <th>#</th>
  <th>Gene ID</th>
  <th>Enzyme Name</th>
  <th>Chromosome</th>
  <th>EC Number</th>
</tr>
      </thead>
      <tbody>
  <?php foreach ($result as $row) { ?>
      <tr>
<td><a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>"><?php echo escape($row["gene_id"]); ?></a></td>
<td><?php echo escape($row["gene_name_short"]); ?></td>
<td><?php echo escape($row["gene_name_full"]); ?></td>
<td><?php echo escape($row["chrom"]); ?></td>
<td><?php echo escape($row["EC_number"]); ?></td>
      </tr>
    <?php } ?>
      </tbody>
  </table>
  <?php } else { ?>
    <p>> No results are available.</p>
  <?php } ?>

<?php
if (isset($_POST['submit'])){?>
  <a href="gene.php">Reset</a>
<?php } ?>
  