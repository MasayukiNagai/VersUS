<!DOCTYPE html>
<html lang="en" ng-app="VersUS-App">
<head>
  <!-- Required meta tags -->
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="An interface for exploring variants of uncertain significance">
  <meta name="author" content="Masayuki Nagai">

  <link rel="stylesheet" href="css/bootstrap.min.css" >
  <link rel="stylesheet" href="css/font-awesome-4.7.0/css/font-awesome.min.css">

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
    a{
      color: #337ab7;
    }
    a[ng-click] {
      cursor: pointer;
    }
    td > a {
      text-decoration: none;
    }
    td > a:hover {
      text-decoration: underline;
    }
    th > a.sortable {
      text-decoration: none;
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

    .alert {
      padding: 15px;
      margin-bottom: 20px;
      border: 1px solid transparent;
      border-radius: 4px;
    }

    .alert-info {
      color: #31708f;
      background-color: #d9edf7;
      border-color: #bce8f1;
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
  <script src="https://ajax.googleapis.com/ajax/libs/angularjs/1.4.9/angular.min.js"></script>
</head>

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

  if(!isset($_GET['sort'])){
    $order = "ORDER BY CADD_score DESC";
  }else{
    $sort_by = $_GET['sort'];
    $desc = $_GET['desc'];
    if($sort_by == 'variation'){
      $order = "ORDER BY ref_pos_alt ";
    } elseif($sort_by == 'cadd'){
      $order = "ORDER BY CADD_score ";
    } elseif($sort_by == 'gnomADAF'){
      $order = "ORDER BY gnomAD_AF ";
    } else{
      $order = "";
    }
    if($desc == 'true'){
      $order .= 'DESC';
    } else{
      $order .= 'ASC';
    }
  }

  $sql2 = "SELECT * FROM Mutation WHERE gene_id = $gene_id
           {$order}
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

<body ng-controller="myCtrl">
<?php require "templates/header2.php"; ?>

<div class="contents">

  <div id="header_mutation">
    <div id="header_mutation_text">
      <h1 class="gene_symbol"><?php echo escape($gene_symbol); ?></h1>
      <h2 class="gene_full_name"><?php echo escape($gene_full_name); ?></h2>
    </div>
  </div>

  <?php require "templates/result_header2.php" ?>
  <hr>

  <?php
  $counter = ($current_page-1) * $num_per_page;
  if ($results && $statement2->rowCount() > 0) { ?>
  <table class="table table-striped table-hover">
    <thead>
      <tr>
        <th><input type="checkbox"></th>
        <th class="count">#</th>
        <th class="variation"><a ng-click="sortType = 'variation'; reverse(); sort()" class="sortable">
            Missense Variation
            <span ng-show="sortType == 'variation' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'variation' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="clinvar_link">Clinvar Link</th>
        <th class="cadd_score"><a ng-click="sortType = 'cadd'; reverse(); sort()" class="sortable">
            CADD Score
            <span ng-show="sortType == 'cadd' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'cadd' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="gnomAD_AF"><a ng-click="sortType = 'gnomADAF'; reverse(); sort()" class="sortable">
            gnomAD AF
            <span ng-show="sortType == 'gnomADAF' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'gnomADAF' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="pdb">PDB</th>
      </tr>
    </thead>
    <tbody>
      <?php foreach ($results as $row) {
        $counter += 1;?>
      <tr>
        <td><input type="checkbox"></td>
        <td class="count"><?php echo escape($counter) ?></td>
        <td class="variation"><?php echo escape($row["ref_pos_alt"]); ?></td>
        <td class="clinvar_link"><a href=<?php echo $row["clinvar_link"] ?>>link</a></td>
        <td class="cadd_score"><?php echo escape(number_format((float)$row["CADD_score"], 1, '.', '')); ?></td>
        <td class="gnomAD_AF"><?php echo escape($row["gnomAD_AF"]); ?></td>
        <td class="pdb"><?php echo escape($row["pdb"]); ?></td>
      </tr>
      <?php } ?>
    </tbody>
  </table>

  <?php } else { ?>
    <p>> No resultss are available.</p>
    <p>> Query: <?php echo escape($sql2) ?></p>
  <?php } ?>

  <script type="text/javascript">

    // sessionStorage.setItem("reverse", true);
    var app = angular.module("VersUS-App", []);
    app.controller("myCtrl", function($scope){
      $scope.sortType = '<?=$_GET['sort']?>';
      $scope.sortReverse = JSON.parse(sessionStorage.getItem('reverse'));
      $scope.sort = function () {tableSort($scope.sortType, JSON.parse(sessionStorage.getItem('reverse')))};
      $scope.reverse = function (){sortReverse()};

      function tableSort(sortType, sortReverse){
        var url = window.location.href;
        var url_query = url.split('?');
        var newurl = ""
        if (url_query.length == 1){
            newurl += url + "?";
        }
        else{
            newurl += url_query[0] + "?";
            var querystring = url_query[url_query.length - 1];
            var queries = querystring.split("&");
            var newQueries = [];
            for (var query of queries){
              if (!query.includes("sort") && !query.includes("desc")){
                  newurl += query + "&";
              }
            }
        }
        newurl += "sort=" + sortType;
        if (sortReverse){
          newurl += "&desc=true";
        }
        location.href = newurl;
      };

      function sortReverse(){
        sessionStorage.reverse = !(JSON.parse(sessionStorage.getItem('reverse')));
      }

    });

  </script>

</div>

<?php require "templates/footer2.php"; ?>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"
    integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
</body>

</html>
