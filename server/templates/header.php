<nav class="navbar navbar-default navbar-expand-md navbar-light bg-light navbar-fixed-top">
  <div class="container">
    <a class="navbar-brand" href="gene.php">VersUS</a>
    <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarSupportedContent"
      aria-controls="navbarSupportedContent" aria-expanded="false" aria-label="Toggle navigation">
      <span class="navbar-toggler-icon"></span>
    </button>

    <div class="collapse navbar-collapse" id="navbarSupportedContent">
      <ul class="navbar-nav me-auto">
        <li class="nav-item active">
          <a class="nav-link" href="gene.php">Home</a>
        </li>
        <li class="nav-item active">
          <a class="nav-link" href="#">About</a>
        </li>
        <li class="nav-item">
          <a class="nav-link" href="tree.html">Tree view</a>
        </li>
      </ul>
      <div class="d-flex justify-content-end align-items-center">
        <form class="d-flex justify-content-end" action="gene.php" method="get" id="search-bar">
          <div class="col-auto">
            <?php
            $search_by = (isset($_GET['search_by'])) ? htmlentities($_GET['search_by']) : '';
            ?>
            <select class="form-select form-select-sm mt-1" id="search_by" name="search_by">
              <option <?php if ($_GET['search_by'] == 'gene') { ?>selected="true" <?php }; ?> value="gene">Gene ID</option>
              <option <?php if ($_GET['search_by'] == 'uniprotID') { ?>selected="true" <?php }; ?> value="uniprotID">Uniprot ID</option>
              <option <?php if ($_GET['search_by'] == 'keywords') { ?>selected="true" <?php }; ?> value="keywords">Keyword</option>
            </select>
          </div>
          <?php
          $search = (isset($_GET['term'])) ? htmlentities($_GET['term']) : '';
          ?>
          <input class="form-control form-control-sm mt-1 mb-1" type="search" name="term" placeholder="Search" value="<?= $search ?>" aria-label="Search">
          <button class="btn btn-primary btn-sm mt-1 mb-1" type="submit">Search</button>
          <!-- <button class="btn btn-default" type="submit">
            <i class="glyphicon glyphicon-search"></i>
          </button> -->
        </form>
        <button type='button' class="btn btn-sm ms-2" ng-class="savedResults.length == 0 ? 'disabled btn-default' : 'btn-primary'" ng-click="checkoutButton(getFileName())">
          <i class="fa fa-archive" aria-hidden="true"></i> <span ng-bind="savedResults.length"></span></button>
        <button type='button' class="btn btn-sm ms-1 btn-default" ng-class="{disabled: savedResults.length == 0}" ng-click="clearSaved()">
          <i class="fa fa-trash" aria-hidden="true"></i> </button>
      </div>
    </div>
  </div>
</nav>
