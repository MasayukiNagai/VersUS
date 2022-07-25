<!DOCTYPE html>
<html lang="en" ng-app="VersUS-App">

<head>
  <!-- Required meta tags -->
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="An interface for exploring variants of uncertain significance">
  <meta name="author" content="Masayuki Nagai">

  <title>VersUS</title>

  <link rel="stylesheet" href="css/bootstrap.min.css" >
  <link rel="stylesheet" href="css/font-awesome-4.7.0/css/font-awesome.min.css">
  <link rel="stylesheet" href="css/style.css">

  <script src="https://ajax.googleapis.com/ajax/libs/angularjs/1.4.9/angular.min.js"></script>
</head>

<?php
require "config.php";
require "common.php";
try{
    $connection = new PDO($dsn, $username, $password, $options);
    $current_page = GetPost('page', 1);
    $num_per_page = 50;
    $start = ($current_page - 1) * $num_per_page;
    $ec = $_GET['ec'];
    $ec_numbers = explode(".", $ec);
    if(count($ec_numbers) == 1){
      $condition = "WHERE ec_1 = {$ec_numbers[0]}";
    }
    elseif(count($ec_numbers) == 2){
      $condition = "WHERE ec_1 = $ec_numbers[0]
                    AND ec_2 = {$ec_numbers[1]}";
    }
    elseif(count($ec_numbers) == 3){
      $condition = "WHERE ec_1 = $ec_numbers[0]
                    AND ec_2 = {$ec_numbers[1]}
                    AND ec_3 = {$ec_numbers[2]}";
    }
    elseif(count($ec_numbers) == 4){
      $condition = "WHERE ec_1 = {$ec_numbers[0]}
                    AND ec_2 = {$ec_numbers[1]}
                    AND ec_3 = {$ec_numbers[2]}
                    AND ec_4 = {$ec_numbers[3]}";
    }
    else{
      echo "Error: " . $ec_numbers;
    }
    if(isset($_GET['sort'])){
      $sort_by = $_GET['sort'];
      $desc = $_GET['desc'];
      if($sort_by == 'gene'){
        $order = "ORDER BY gene_symbol ";
      } elseif($sort_by == 'uniprot'){
        $order = "ORDER BY uniprot_id ";
      } elseif($sort_by == 'name'){
        $order = "ORDER BY gene_full_name ";
      }elseif($sort_by == 'vus'){
        $order = "ORDER BY num_vus ";
      } elseif($sort_by == 'cadd'){
        $order = "ORDER BY max_cadd ";
      } elseif($sort_by == 'ec'){
        $order = "ORDER BY EC_number ";
      } else{
        $order = "";
      }
      if($desc == 'true'){
        $order .= 'DESC';
      } else{
        $order .= 'ASC';
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

<body ng-controller="myCtrl" ng-init="init()">
<?php require "templates/header.php"; ?>

<div class="contents">

  <?php require "templates/pagination.php" ?>

  <hr>

  <div class="table">
  <?php
  $counter = ($current_page-1) * $num_per_page;
  if ($results && $statement->rowCount() > 0) { ?>
    <table class="table table-striped table-hover">
      <thead>
        <tr>
        <th class="count">#</th>
        <th class="gene_id"><a ng-click="sortType = 'gene'; reverse(); sort()" class="sortable">
            Gene
            <span ng-show="sortType == 'gene' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'gene' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="enzyme_name"><a ng-click="sortType = 'name'; reverse(); sort()" class="sortable">
            Enzyme Name
            <span ng-show="sortType == 'name' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'name' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="uniprot_id"><a ng-click="sortType = 'uniprot'; reverse(); sort()" class="sortable">
            Uniprot ID
            <span ng-show="sortType == 'uniprot' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'uniprot' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="num_vus"><a ng-click="sortType = 'vus'; reverse(); sort()" class="sortable">
            # Missense VUS
            <span ng-show="sortType == 'vus' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'vus' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="cadd_score"><a ng-click="sortType = 'cadd'; reverse(); sort()" class="sortable">
            Highest CADD Score
            <span ng-show="sortType == 'cadd' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'cadd' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="EC_number"><a ng-click="sortType = 'ec'; reverse(); sort()" class="sortable">
            EC #
            <span ng-show="sortType == 'ec' && sortReverse" class="fa fa-caret-down"></span>
            <span ng-show="sortType == 'ec' && !sortReverse" class="fa fa-caret-up"></span>
          </a></th>
        <th class="alphafold">AlphaFold Protein Structure</th>
        </tr>
      </thead>
      <tbody>
      <?php foreach ($results as $row) {
        $counter += 1; ?>
        <tr>
          <td class="count"><?php echo escape($counter) ?></td>
          <td class="gene_id"><a href="variant.php?gene_id=<?php echo $row["gene_id"] ?>"><?php echo escape($row["gene_symbol"]); ?></a></td>
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
    <p>> Query: <?php echo escape($sql) ?></p>
  <?php } ?>
  </div>

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

      // $scope.passMutationInfo = function (){
      //   var php_results = <?= json_encode($results); ?>;
      //   for (var i = 0; i < php_results.length; i++){
      //     value = php_results[i];
      //     $scope.results.push({
      //       'mutation_id': value['mutation_id'],
      //       'gene_id': value['gene_id'],
      //       'ref': value['ref'],
      //       'pos': value['pos'],
      //       'alt': value['alt'],
      //       'accession': value['accession'],
      //       'CADD_score': value['CADD_score'],
      //       'gnomAD_AF': value['gnomAD_AF'],
      //       'pdb': value['pdb'],
      //       'fasta_id': value['fasta_id']
      //     })
      //   }
      // };

      // $scope.passFastaInfo = function () {
      //   var php_fasta = <?= json_encode($fasta_seqs); ?>;
      //   for (var i = 0; i < php_fasta.length; i++){
      //     value = php_fasta[i];
      //     $scope.fastaSeqs.push({
      //       'fasta_id': value['fasta_id'],
      //       'NP_accession': value['NP_accession'],
      //       'fasta': value['fasta']
      //     })
      //   }
      // };

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
        $scope.initSortReverse();
        $scope.initSavedResults();
      };

    });
  </script>

</div>

<?php require "templates/footer.php"; ?>
  <!-- <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
  <script src="js/bootstrap.min.js"></script>
  <link href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" rel="stylesheet" type="text/css" />
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script> -->
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"
    integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
</body>

</html>
