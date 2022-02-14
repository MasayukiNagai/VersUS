<?php
# If empty term is entered, redirect to the main page
if(isset($_GET['term']) && (empty($_GET['term']) || ctype_space($_GET['term']))){
    header('Location: gene.php');
    exit();
}
?>

<?php require "templates/header.php"; ?>
<div class="content">

<?php
require "config.php";
require "common.php";
try{
    $connection = new PDO($dsn, $username, $password, $options);
    $current_page = GetPost('page', 1);
    $num_per_page = 50;
    $start = ($current_page - 1) * $num_per_page;
    if(!isset($_GET['term'])){
        $condition = "";
        $order = "ORDER BY gene_symbol ASC";
    }
    else{
        $search_by = $_GET['search_by'];
        $term = rtrim($_GET['term']);
        if ($search_by == 'gene'){
          // search by Gene ID
          $condition = "WHERE gene_symbol LIKE :keyword";
          $order = "ORDER BY gene_symbol ASC";
          $keyword = "$term%";
        }
        elseif($search_by == 'uniprotID'){
          // search by Uniprot ID
          $condition = "WHERE uniprot_id = :keyword";
          $order = "ORDER BY uniprot_id ASC";
          $keyword = "$term";
        }
        else{
          // search by keywords
          $condition = "WHERE gene_full_name LIKE :keyword
                           OR EC_number LIKE :keyword";
          $order = "";
          $keyword = "%$term%";
        }
    }
    $limit = "LIMIT :start, :num_per_page";
    # query to get gene items
    $sql = get_query($condition, $order, $limit);
    $statement = $connection->prepare($sql);
    $statement->bindParam(':start', $start, PDO::PARAM_INT);
    $statement->bindParam(':num_per_page', $num_per_page, PDO::PARAM_INT);
    if(isset($_GET['term'])){
        $statement->bindParam(':keyword', $keyword, PDO::PARAM_STR);
    }
    $statement->execute();
    $results = $statement->fetchAll();

    # query to get the number of results
    $sql2 = "SELECT COUNT(*) FROM Gene {$condition}";
    $statement2 = $connection->prepare($sql2);
    if(isset($_GET['term'])){
        $statement2->bindParam(':keyword', $keyword, PDO::PARAM_STR);
    }
    $statement2->execute();
    $num_results = $statement2->fetch()[0];
    $total_page = ceil($num_results/$num_per_page);
}catch(PDOException $error) {
    echo $error->getMessage();
}
?>


<?php require "templates/result_header.php" ?>

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
        <th class="alphafold">AlphaFold Protein Structure</th>
        </tr>
      </thead>
      <tbody>
      <?php foreach ($results as $row) {
        $counter += 1; ?>
        <tr>
        <td class="count"><?php echo escape($counter) ?></td>
        <td class="gene_id"><a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>"><?php echo escape($row["gene_symbol"]); ?></a></td>
        <td class="enzyme_name"><?php echo escape($row["gene_full_name"]); ?></td>
        <td class="uniprot_id"><a href=<?php echo get_uniprot_url($row["uniprot_id"]) ?>><?php echo escape($row["uniprot_id"]) ?></a></td>
        <td class="num_vus"><?php echo escape($row["num_vus"]); ?></td>
        <td class="cadd_score"><?php echo escape(number_format((float)$row["max_cadd"], 1, '.', '')); ?></td>
        <td class="EC_number"><?php echo escape($row["EC_number"]); ?></td>
        <td class="alphafold"><a href=<?php echo get_alphafold_url($row["uniprot_id"]) ?>>link</a></td>
        </tr>
      <?php } ?>
      </tbody>
    </table>

<?php } else { ?>
  <p>> No results are available.</p>
<?php } ?>

</div>
<?php require "templates/footer.php"; ?>
