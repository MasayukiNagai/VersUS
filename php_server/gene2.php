<!DOCTYPE html>
<html lang="en">

<head>
  <!-- Required meta tags -->
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="An interface for exploring variants of uncertain significance">
  <meta name="author" content="Masayuki Nagai">

  <link href="css/bootstrap.min.css" rel="stylesheet">

  <title>VersUS(beta)</title>

  <style type="text/css">

    .contents{
      font-size: 15px;
      width: 90%;
      /* margin-left: auto;
      margin-right: auto; */
      margin: 0px auto 50px;
    }
    hr{
      color: grey;
    }
    /*
    Result Header
    */
    .contents_header{
        margin-top: 5px;
        margin-bottom: 30px;
        display: inline-block;
        /* height: 150px; */
        /* border: solid 1px black; */
    }

    .btn-group-xs > .btn, .btn-xs {
    padding: 1px 5px;
    font-size: 12px;
    line-height: 1.5;
    border-radius: 3px;
    }

    .num_items{
        position: absolute;
        left: 5%;
        display: inline-block;
        margin: 0 auto;
    }

    .pageforms{
        position: absolute;
        right: 5%;
        display: inline-block;
        vertical-align: bottom;
        margin: 5px auto;
    }

    .pageform{
        display: inline-block;
    }

    .page_button, #jump_button{
        cursor: pointer;
    }

    .page_button:disabled{
        cursor: auto;
    }

    #page_number{
      /* height:15px; */
      width: 50px;
      padding:0 auto;
      position:relative;
      margin: 0 2px;
      /* top:0;  */
      border-radius:2px;
      border: 1px solid grey;
      /* outline:0; */
      background:white;
      border: solid 1px black;
    }

    footer{
      width: 90%;
      margin: 50px auto 50px;
      color: grey;
      font-size: 14px;
    }

  </style>

  <!-- <style type="text/css">
    .navbar-custom {
    background-color: #633974 ;
    }
    /* change the brand and text color */
    .navbar-custom .navbar-brand,
    .navbar-custom .navbar-text {
    color: rgba(255,255,255,.8);
    }

    /* change the color of active or hovered links */
    .navbar-custom .nav-item.active .nav-link,
    .navbar-custom .nav-item:hover .nav-link {
    color: #ffffff;
    }
  </style> -->
</head>

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

<body>
<?php require "templates/header2.php"; ?>

<div class="contents">

  <?php require "templates/result_header2.php" ?>

  <hr>

  <div class="table">
  <?php
  $counter = ($current_page-1) * $num_per_page;
  if ($results && $statement->rowCount() > 0) { ?>
    <table class="table table-striped table-hover">
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
          <td class="gene_id"><a href="mutation2.php?gene_id=<?php echo $row["gene_id"] ?>"><?php echo escape($row["gene_symbol"]); ?></a></td>
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
</div>

<?php require "templates/footer2.php"; ?>
  <!-- <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
  <script src="js/bootstrap.min.js"></script>
  <link href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" rel="stylesheet" type="text/css" />
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script> -->
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"
    integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
</body>

</html>
