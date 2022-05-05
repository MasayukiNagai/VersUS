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
    .btn-primary, .btn-primary.disabled {
      color: #fff;
      background-color: #337ab7;
      border-color: #2e6da4;
    }
    .btn-primary:hover {
      color: #fff;
      background-color: #286090;
      border-color: #204d74;
    }
    .btn-primary.focus, .btn-primary:focus {
      color: #fff;
      background-color: #286090;
      border-color: #122b40;
    }
    .btn-default {
      color: #333;
      background-color: #fff;
      border-color: #ccc;
    }
    .btn.disabled{
      cursor: not-allowed;
      filter: alpha(opacity=65);
      -webkit-box-shadow: none;
      box-shadow: none;
      opacity: .65;
    }
    /*
    Result Header
    */
    .contents_header{
        margin-top: 5px;
        margin-bottom: 5px;
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
        /* position: absolute; */
        left: 5%;
        display: inline-block;
        margin: 0 auto;
    }

    .pageforms{
        /* position: static; */
        left: 5%;
        /* display: inline-block; */
        margin: 3px auto;
        /* vertical-align: bottom; */
        /* border: solid 1px black; */
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
      $order = "ORDER BY pos ";
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

  $sql4 = "SELECT fasta_id, NP_accession, fasta FROM Fasta RIGHT JOIN
          (SELECT DISTINCT fasta_id FROM Mutation
           WHERE gene_id = $gene_id) AS m USING(fasta_id)";
  $statement4 = $connection->prepare($sql4);
  $statement4->execute();
  $fasta_seqs = $statement4->fetchAll();

}catch(PDOException $error) {
  echo $sql2 . "<br>" . $error->getMessage();
}?>

<body ng-controller="myCtrl">
<?php require "templates/header2.php"; ?>

<div class="contents" ng-init="init()">

  <div id="header_mutation">
    <div id="header_mutation_text">
      <h1 class="gene_symbol"><?php echo escape($gene_symbol); ?></h1>
      <h2 class="gene_full_name"><?php echo escape($gene_full_name); ?></h2>
    </div>
  </div>

  <?php require "templates/result_header2.php" ?>
  <hr>
  <button type='button' class="btn btn-primary" ng-class="{disabled: selectedItems == 0}" ng-click="saveDatasets()">
    Add <span ng-bind="selectedItems">0</span> to collection</button>
  <button type='button' class="btn" ng-class="savedResults.length == 0 ? 'disabled btn-default' : 'btn-primary'" ng-click="checkoutButton('test.fasta')">
    <i class="fa fa-archive" aria-hidden="true"></i> Download <span ng-bind="savedResults.length"></span> fasta seqs</button>
  <button type='button' class="btn btn-default" ng-show="savedResults.length > 0" ng-click="clearSaved()">
    Clear the collection </button>
  <hr>
  <?php
  $counter = ($current_page-1) * $num_per_page;
  if ($results && $statement2->rowCount() > 0) { ?>
  <table class="table table-striped table-hover">
    <thead>
      <tr>
        <th><input type="checkbox" ng-model="selectAll" ng-click="checkAll()"></th>
        <!-- <th class="count">#</th> -->
        <th class="variation"><a ng-click="sortType = 'variation'; reverse(); sort()" class="sortable">
            Missense Variation
            <span ng-show="sortType == 'variation' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'variation' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="clinvar_link">Clinvar Accession</th>
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
      <tr ng-repeat="x in results" ng-click="rowClicked(x)">
        <td><input type="checkbox" ng-checked="x.selected" ng-click="toggleRow($event, x)"></td>
        <td>{{ x.ref + x.pos + x.alt}}</td>
        <td><a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{{ x.accession }}"> {{ x.accession }}</a></td>
        <td>{{ x.CADD_score }}</td>
        <td>{{ x.gnomAD_AF }}</td>
        <td>{{ x.pdb }}</td>
      </tbody>
  </table>

  <?php } else { ?>
    <p>> No resultss are available.</p>
    <p>> Query: <?php echo escape($sql2) ?></p>
  <?php } ?>

  <script type="text/javascript">

    // sessionStorage.setItem("reverse", true);
    var app = angular.module("VersUS-App", []);
    var aaMapThreeToOne = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                           'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                           'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                           'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}
    app.controller("myCtrl", function($scope){
      $scope.gene_symbol = <?= json_encode($gene_symbol); ?>;
      $scope.sortType = '<?= $_GET['sort'] ?>';
      // $scope.sortReverse = JSON.parse(sessionStorage.reverse);
      $scope.sort = function () {tableSort($scope.sortType, JSON.parse(sessionStorage.reverse))};
      $scope.reverse = function (){sortReverse()};

      $scope.results = [];
      // $scope.savedResults = [];
      $scope.fastaSeqs = []

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
        sessionStorage.reverse = !(JSON.parse(sessionStorage.reverse));
      };

      $scope.passMutationInfo = function (){
        var php_results = <?= json_encode($results); ?>;
        for (var i = 0; i < php_results.length; i++){
          value = php_results[i];
          $scope.results.push({
            'mutation_id': value['mutation_id'],
            'gene_id': value['gene_id'],
            'ref': value['ref'],
            'pos': value['pos'],
            'alt': value['alt'],
            'accession': value['accession'],
            'CADD_score': value['CADD_score'],
            'gnomAD_AF': value['gnomAD_AF'],
            'pdb': value['pdb'],
            'fasta_id': value['fasta_id']
          })
        }
      };

      $scope.passFastaInfo = function () {
        var php_fasta = <?= json_encode($fasta_seqs); ?>;
        for (var i = 0; i < php_fasta.length; i++){
          value = php_fasta[i];
          $scope.fastaSeqs.push({
            'fasta_id': value['fasta_id'],
            'NP_accession': value['NP_accession'],
            'fasta': value['fasta']
          })
        }
      };

      $scope.toggleRow = function ($event, obj) {
        $event.stopPropagation();
        obj.selected = !obj.selected;
      };

      $scope.rowClicked = function (obj) {
        obj.selected = !obj.selected;
      };

      $scope.checkAll = function () {
        angular.forEach($scope.results, function (item) {
          item.selected = $scope.selectAll;
        });
      };

      $scope.$watch('results', function (items) {
        var selectedItems = 0;
        angular.forEach(items, function (item) {
          selectedItems += item.selected ? 1 : 0;
        });
        $scope.selectedItems = selectedItems;
      }, true);

      $scope.saveDatasets = function () {
        angular.forEach($scope.results, function (item) {
          if (item.selected) {
            var alreadyAdded = false;
            for (var i = 0; i < $scope.savedResults.length; i++) {
              if (item['mutation_id'] == $scope.savedResults[i]['mutation_id']) {
                alreadyAdded = true;
              }
            }
            if (!alreadyAdded) {
              // $scope.savedResults.push(angular.copy(item));
              // fasta seq instead of fasta_id?
              // gene name instead of gene_id?
              for(var i = 0; i < $scope.fastaSeqs.length; i++) {
                if (item['fasta_id'] == $scope.fastaSeqs[i]['fasta_id']) {
                  var index = i;
                  break;
                }
              }
              $scope.savedResults.push({
                'mutation_id': item['mutation_id'],
                'fasta_id': item['fasta_id'],
                'ref': item['ref'],
                'pos': item['pos'],
                'alt': item['alt'],
                'gene_id': item['gene_id'],
                'gene_symbol': $scope.gene_symbol,
                'NP_accession': $scope.fastaSeqs[index]['NP_accession'],
                'fasta': $scope.fastaSeqs[index]['fasta']
              })
            }
          }
          sessionStorage.savedResults = JSON.stringify($scope.savedResults);
        });
      };

      $scope.clearSaved = function () {
        $scope.savedResults = [];
        sessionStorage.removeItem("savedResults");
      }

      $scope.checkoutButton = function (filename) {
        var fasta = $scope.prepFasta();
        $scope.downloadFile(fasta, filename);
        // $scope.savedResults = [];
      };

      $scope.prepFasta = function () {
        var contents = [];
        for (var i = 0; i < $scope.savedResults.length; i++){
          var line = [];
          line.push('>' + $scope.savedResults[i]['NP_accession']);
          line.push($scope.savedResults[i]['ref']+$scope.savedResults[i]['pos']+$scope.savedResults[i]['alt']);
          contents.push(line.join('\t'));
          fasta = $scope.savedResults[i]['fasta'];
          pos = $scope.savedResults[i]['pos'];
          ref = aaMapThreeToOne[$scope.savedResults[i]['ref']];
          alt = aaMapThreeToOne[$scope.savedResults[i]['alt']];
          mtFasta = fasta.modifyFasta(Number(pos), ref, alt);
          contents.push(mtFasta + '\n');
        }
        return contents.join('\n');
      }

      $scope.downloadFile = function (content, filename) {
        var blob = new Blob([content], { type: "text/plain;charset=utf-8" });
        // saveAs(blob, filename);
        var downloadLink = document.createElement("a");
        downloadLink.download = filename;
        downloadLink.href = window.URL.createObjectURL(blob);
        downloadLink.style.display = "none";
        document.body.appendChild(downloadLink);
        downloadLink.click();
      };

      String.prototype.modifyFasta = function(pos, ref, alt) {
        if (fasta[pos-1] == ref) {
          return this.substring(0, pos-1) + alt + this.substring(pos-1 + alt.length);
        }
        else{
          console.log(ref + pos + alt);
          return null;
        }
      };

      $scope.initSortReverse = function () {
        if (sessionStorage.reverse == null){
          sessionStorage.reverse = false;
        }
        else{
          $scope.sortReverse = JSON.parse(sessionStorage.reverse);
        }
      }

      $scope.initSavedResults = function () {
        if (sessionStorage.savedResults == null) {
          // sessionStorage.savedResults = [];
          $scope.savedResults = [];
        }
        else {
          // console.log(sessionStorage.savedResults);
          $scope.savedResults = JSON.parse(sessionStorage.savedResults);
        }
      };

      $scope.init = function(){
        console.log('Call init function');
        $scope.passMutationInfo();
        $scope.passFastaInfo();
        $scope.initSortReverse();
        $scope.initSavedResults();
      };

    });

  </script>

</div>

<?php require "templates/footer2.php"; ?>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"
    integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
</body>

</html>
