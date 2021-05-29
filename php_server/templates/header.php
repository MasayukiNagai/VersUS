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
      <!-- <div id="headinfo">
        <a href="gene.php?page=1">Home</a>
      </div> -->
      <div id="header">
        <h1 id="title-image"><a href="gene.php?page=1"><img src="images/versus_logo.png" alt="versus" /></a></h1>
        <div id="queryContainer">
          <form method="post" id="search-bar">
            <select id="search_by" name="search_by">
              <option value="gene_name_short">Gene ID</option>
              <option value="uniprot_id">Uniprot ID</option>
              <option value="keywords">Keyword</option>
            </select>
            <input type="text" id="textbox" name="keyword">
            <input type="submit" name="submit" id="sbtn" value="Search">
          </form>
        </div>
      </div>
    </div>
