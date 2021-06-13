<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <meta http-equiv="x-ua-compatible" content="ie=edge" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>VersUS</title>
    <link rel="stylesheet" href="css/style.css" />
  </head>

  <body>
    <div>
      <div id="header">
        <h1 id="title-image"><a href="gene.php?page=1"><img src="images/versus_logo.png" alt="versus" /></a></h1>
        <div id="queryContainer">
          <form action="gene.php", method="get" id="search-bar">
            <select id="search_by" name="search_by">
              <option value="gene">Gene ID</option>
              <option value="uniprotID">Uniprot ID</option>
              <option value="keywords">Keyword</option>
            </select>
            <?php 
            $search = (isset($_GET['term'])) ? htmlentities($_GET['term']) : '';
            ?>
            <input type="text" id="textbox" name="term" value="<?= $search ?>">
            <!-- <input type="submit" name="submit" id="sbtn" value="Search"> -->
            <button id="sbtn">Search</button>
          </form>
        </div>
      </div>
    </div>
