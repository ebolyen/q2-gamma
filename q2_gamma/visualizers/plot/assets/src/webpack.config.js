const path = require('path');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const HtmlWebpackIncludeAssetsPlugin = require('html-webpack-include-assets-plugin');

module.exports = {
  entry: './app/index.js',
  output: {
    filename: 'bundle.js',
    path: path.resolve(__dirname, '../dist')
  },
  mode: 'development',
  module: {
    rules: [
      {
        test: /\.rs$/,
        use: [{
          loader: 'wasm-loader'
        }, {
          loader: 'rust-native-wasm-loader',
          options: {
            release: true
          }
        }]
      }
    ]
  },
  plugins: [
    new HtmlWebpackPlugin({
      title: 'gamma plot'
    }),
    new HtmlWebpackIncludeAssetsPlugin({
      assets: ['data/packed_table.jsonp', 'data/tree.jsonp'],
      append: true,
      jsExtensions: ['.jsonp', '.js']
    })
  ]
}

