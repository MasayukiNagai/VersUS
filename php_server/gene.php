<?php 
# If empty term is entered, redirect to the main page
if(isset($_GET['term']) && (empty($_GET['term']) || ctype_space($_GET['term']))){
    header('Location: gene.php');
    exit();
}
?>

<?php require "templates/header2.php"; ?>
<div class="content">

<?php
require "config.php";
require "common.php";
try{
    $connection = new PDO($dsn, $username, $password, $options);
    $current_page_count = Get('page', 1);
    $num_per_page = 50;
    $start = ($current_page_count - 1) * $num_per_page;
    if(!isset($_GET['term'])){
        $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                    COUNT(m.mutation_id) AS num_vus, 
                    MAX(m.CADD_score) AS max_cadd, 
                    g.EC_number
                FROM (SELECT * FROM Gene 
                      ORDER BY gene_name_short ASC LIMIT :start, $num_per_page) as g
                LEFT JOIN Mutation AS m USING(gene_id)
                GROUP BY g.gene_id";
        $statement = $connection->prepare($sql);
        $statement->bindParam(':start', $start, PDO::PARAM_INT);
        $statement->execute();
        $results = $statement->fetchAll();

        $sql2 = "SELECT COUNT(*) FROM Gene";
        $statement2 = $connection->prepare($sql2);
        $statement2->execute();
        $result2 = $statement2->fetch();
        $total_page_count = ceil($result2[0]/$num_per_page);        
    }
    else{
        $search_by = $_GET['search_by'];
        $term = rtrim($_GET['term']);
        if ($search_by == 'gene'){
          // search by Gene ID
          $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                         COUNT(m.mutation_id) AS num_vus, 
                         MAX(m.CADD_score) AS max_cadd, 
                         g.EC_number
                  FROM (SELECT * FROM Gene WHERE gene_name_short LIKE :keyword) AS g
                  LEFT JOIN Mutation AS m USING(gene_id)
                  GROUP BY g.gene_id
                  LIMIT 50";
          $sql2 = "SELECT COUNT(*) FROM Gene WHERE gene_name_short LIKE :keyword";
          $keyword = "$term%";
        }
        elseif($search_by == 'uniprotID'){
          // search by Uniprot ID
          $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                         COUNT(m.mutation_id) AS num_vus, 
                         MAX(m.CADD_score) AS max_cadd, 
                         g.EC_number
                  FROM (SELECT * FROM Gene WHERE uniprot_id = :keyword) AS g
                  LEFT JOIN Mutation AS m USING(gene_id)
                  GROUP BY g.gene_id
                  LIMIT 50";
          $sql2 = "SELECT COUNT(*) FROM Gene WHERE uniprot_id = :keyword";
          $keyword = "$term";
        }
        else{ 
          // search by keywords
          $sql = "SELECT g.gene_id, g.gene_name_short, g.gene_name_full, 
                         COUNT(m.mutation_id) AS num_vus, 
                         MAX(m.CADD_score) AS max_cadd, 
                         g.EC_number
                  FROM (SELECT * FROM Gene
                        WHERE gene_name_short LIKE :keyword
                           OR gene_name_full LIKE :keyword
                           OR EC_number LIKE :keyword) AS g
                  LEFT JOIN Mutation AS m USING(gene_id)
                  GROUP BY g.gene_id
                  LIMIT 50";
          $sql2 = "SELECT COUNT(*) FROM Gene WHERE gene_name_short LIKE :keyword
                                             OR gene_name_full LIKE :keyword
                                             OR EC_number LIKE :keyword";
          $keyword = "%$term%";
        }
        $statement = $connection->prepare($sql);
        $statement->bindParam(':keyword', $keyword, PDO::PARAM_STR);
        $statement->execute();
        $results = $statement->fetchAll();
        
        $statement2 = $connection->prepare($sql2);
        $statement2->bindParam(':keyword', $keyword, PDO::PARAM_STR);
        $statement2->execute();
        $value2 = $statement2->fetch();
        $total_page_count = ceil($value2[0]/$num_per_page);
    }
}catch(PDOException $error) {
    echo $error->getMessage();
}
?>


<?php
if ($results && $statement->rowCount() > 0) { ?>
    <table>
    　<thead>
        <tr>
        <th class="gene_id">Gene ID</th>
        <th class="enzyme_name">Enzyme Name</th>
        <th class="num_vus"># Missense VUS</th>
        <th class="cadd_score">Highest CADD score</th>
        <th class="EC_number">EC #</th>
        </tr>
      </thead>
      <tbody>
      <?php foreach ($results as $row) { ?>
        <tr>
        <td class="gene_id"><a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>&page=1"><?php echo escape($row["gene_name_short"]); ?></a></td>
        <td class="enzyme_name"><?php echo escape($row["gene_name_full"]); ?></td>
        <td class="num_vus"><?php echo escape($row["num_vus"]); ?></td>
        <td class="cadd_score"><?php echo escape($row["max_cadd"]); ?></td>
        <td class="EC_number"><?php echo escape($row["EC_number"]); ?></td>
        </tr>
      <?php } ?>
      </tbody>
    </table>

<?php require "templates/pagination.php";?>


<?php } else { ?>
  <p>> No results are available.</p>
<?php } ?>

</div>
<?php require "templates/footer.php"; ?>