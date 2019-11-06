


<!DOCTYPE HTML>
<html>

<head>
    <meta charset="utf-8">

    <title>JupyterHub</title>
    <meta http-equiv="X-UA-Compatible" content="chrome=1">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    
    <link rel="stylesheet" href="/hub/static/css/style.min.css?v=883a1cf654bd9f5cdc5ca7c91efe7d71" type="text/css"/>
    
    <script src="/hub/static/components/requirejs/require.js?v=f0cc8bbb2fcef87fc194fecbb632fcfa" type="text/javascript" charset="utf-8"></script>
    <script src="/hub/static/components/jquery/dist/jquery.min.js?v=a09e13ee94d51c524b7e2a728c7d4039" type="text/javascript" charset="utf-8"></script>
    <script src="/hub/static/components/bootstrap/dist/js/bootstrap.min.js?v=2f34b630ffe30ba2ff2b91e3f3c322a1" type="text/javascript" charset="utf-8"></script>
    <script>
      require.config({
          
          urlArgs: "v=20191001085403",
          
          baseUrl: '/hub/static/js',
          paths: {
            components: '../components',
            jquery: '../components/jquery/dist/jquery.min',
            bootstrap: '../components/bootstrap/dist/js/bootstrap.min',
            moment: "../components/moment/moment",
          },
          shim: {
            bootstrap: {
              deps: ["jquery"],
              exports: "bootstrap"
            },
          }
      });
    </script>

    <script type="text/javascript">
      window.jhdata = {
        base_url: "/hub/",
        prefix: "/",
        
        
        admin_access: false,
        
        
        options_form: false,
        
      }
    </script>

    
    

</head>

<body>

<noscript>
  <div id='noscript'>
    JupyterHub requires JavaScript.<br>
    Please enable it to proceed.
  </div>
</noscript>


  <nav class="navbar navbar-default">
    <div class="container-fluid">
      <div class="navbar-header">
        <span id="jupyterhub-logo" class="pull-left"><a href="/hub/"><img src='/hub/logo' alt='JupyterHub' class='jpy-logo' title='Home'/></a></span>
        <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#thenavbar" aria-expanded="false">
          <span class="sr-only">Toggle navigation</span>
          <span class="icon-bar"></span>
          <span class="icon-bar"></span>
          <span class="icon-bar"></span>
        </button>
      </div>

      <div class="collapse navbar-collapse" id="thenavbar">
        
        <ul class="nav navbar-nav navbar-right">
          
            <li>
              

            </li>
          
        </ul>
      </div>

      
      
    </div>
  </nav>











<div id="login-main" class="container">

<div class="service-login">
  <a role="button" class='btn btn-jupyter btn-lg' href='/hub/saml2_auth/login?next=%2Fhub%2Fapi%2Foauth2%2Fauthorize%3Fclient_id%3Djupyterhub-user-sch14%26redirect_uri%3D%252Fuser%252Fsch14%252Foauth_callback%26response_type%3Dcode%26state%3DeyJ1dWlkIjogIjE3ZmZkNWVkZWMzNDQ3ZGFiZTRmNTJkN2FiNTgxMWUwIiwgIm5leHRfdXJsIjogIi91c2VyL3NjaDE0L2ZpbGVzL3ZvcGVyYXRlLnB5In0'>
    Sign in with college ID
  </a>
</div>

</div>







<div class="modal fade" id="error-dialog" tabindex="-1" role="dialog" aria-labelledby="error-label" aria-hidden="true">
  <div class="modal-dialog">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal"><span aria-hidden="true">&times;</span><span class="sr-only">Close</span></button>
        <h4 class="modal-title" id="error-label">Error</h4>
      </div>
      <div class="modal-body">
        
  <div class="ajax-error">
    The error
  </div>

      </div>
      <div class="modal-footer">
        <button type="button" class="btn btn-default" data-dismiss="modal">Cancel</button>
        <button type="button" class="btn btn-primary" data-dismiss="modal" data-dismiss="modal">OK</button>
      </div>
    </div>
  </div>
</div>





<script>
if (window.location.protocol === "http:") {
  // unhide http warning
  var warning = document.getElementById('insecure-login-warning');
  warning.className = warning.className.replace(/\bhidden\b/, '');
}
</script>



</body>

</html>