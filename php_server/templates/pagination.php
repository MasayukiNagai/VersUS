<?php
if(!isset($_GET['term'])){?>
    <div class="pagination">
        <?php 
        if($current_page_count > 1){?>
            <a href="gene.php?page=<?php echo $current_page_count-1 ?>">&laquo;</a>
        <?php }
        $page_count = 1;
        while($page_count <= $total_page_count){?>
            <a href="gene.php?page=<?php echo $page_count ?>"><?php echo escape($page_count) ?></a>
            <?php $page_count++;
        } 
        if($current_page_count < $total_page_count){?>
            <a href="gene.php?page=<?php echo $current_page_count+1 ?>">&raquo;</a>
        <?php }?> 
    </div>
<?php }else{ ?>
    <a href="gene.php">Reset</a>
    <div class="pagination">
        <?php 
        if($current_page_count > 1){?>
            <a href="gene.php?search_by=<?php echo htmlentities($search_by) ?>&term=<?php echo htmlentities($term) ?>&page=<?php echo $current_page_count-1 ?>">&laquo;</a>
        <?php }
        $page_count = 1;
        while($page_count <= $total_page_count){?>
            <a href="gene.php?search_by=<?php echo htmlentities($search_by) ?>&term=<?php echo htmlentities($term) ?>&page=<?php echo $page_count ?>"><?php echo escape($page_count) ?></a>
            <?php $page_count++;
        } 
        if($current_page_count < $total_page_count){?>
            <a href="gene.php?search_by=<?php echo htmlentities($search_by) ?>&term=<?php echo htmlentities($term) ?>&page=<?php echo $current_page_count+1 ?>">&raquo;</a>
        <?php }?> 
    </div>
<?php }?>