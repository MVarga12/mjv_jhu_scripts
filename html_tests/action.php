<?php
if (isset($_POST['email'])) {
    $email_to = "matthewjvarga@gmail.com";
    $email_subject = "new document request";

    function died($error) {
        echo "An error occured: <br>";
        echo $error. "<br>";
        die();
    }

    if(!isset($_POST['user_name']) ||
        (!isset($_POST['email']))) {
        died("Please enter your name and email.");
    }

    $error_message = $first = $last = $tit = $journal = "";
    $user_name = htmlspecialchars($_POST['user_name']);
    $email = htmlspecialchars($_POST['email']);
    $email_from = $_POST['email'];
    $email_exp = '/^[A-Za-z0-9._%]+@[A-Za-z0-9]+\.[A-Za-z]{2,4}%/';
    if(!preg_match($email_exp,$email_from)) {
        $error_message .= 'Please enter a valid email';
    }

    $first = htmlspecialchars($_POST['FirstName']);
    $last = htmlspecialchars($_POST['LastName']);
    $tit = htmlspecialchars($_POST['Title']);
    $journal = htmlspecialchars($_POST['Journal']);

    $email_message = "Document request details below.\n\n";
    function clean_string($string) {
        $bad = array("content-type","bcc:","to","cc:","href");
        return str_replace($bad,"",$string);
    }

    $email_message .= "Requester: ".clean_string($user_name)."\n";
    #$email_message .= "Document Type: ".clean_string($Type)."\n";
    #$email_message .= "First Author: ".clean_string($Last).", ".clean_string($First)."\n";
    #$email_message .= "Title: ".clean_string($Title)."\n";
    #$email_message .= "Jouenal: ".clean_string($Journal)."\n";

    $headers = 'From: '.$email_from."\r\n".
        'Reply-To: '.$email_from."\r\n".
        'X-Mailer: PHP/' . phpversion();
    @mail($email_to, $email_subject, $email_message, $headers);
?>
    Thank you for requesting a document. It will be added shortly.

<?php
}
    die();
?>
